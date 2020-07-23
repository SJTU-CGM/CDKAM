/**
# CDKAM: a metagenomic classification tool using discriminative k-mers and approximate matching strategy
# Copyright 2019-2020
# Department of Bioinformatics and Biostatistics, Shanghai Jiao Tong University
# Contact information: buikien.dp@sjtu.edu.cn, ccwei@sjtu.edu.cn
#
# Function: Reading sequences and carrying out classification
*/


#include <omp.h>
#include "seqreader.h"
#define LL long long
#define ULL unsigned long long
#define fi first
#define se second
#define FOR(i,a,b) for(size_t i=a;i<=b;i++)
#define FO(i,a,b) for(size_t i=a;i<b;i++)
#define DEBUG(a) {cerr << #a << ": " << (a) << endl; fflush(stderr); }


template<class T> string i2s(T x) {ostringstream o; o << x; return o.str();}
template<class T> int getbit(T s, int i) { return (s >> i) & 1; }
template<class T> T onbit(T s, int i) { return s | (T(1) << i); }
template<class T> T offbit(T s, int i) { return s & (~(T(1) << i)); }
template<class T> int cntbit(T s) { return __builtin_popcount(s);}
typedef pair<int, int> II;


/***************************************************************/
#include <sys/time.h>
#include <sys/resource.h>
inline void printRam() {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    cerr << "Max ram (in kilobytes): " << ru.ru_maxrss << endl;
}

namespace Time{
    double start_time, time_limit;
    static double last_call = 0;
    int get_time_calls = 0;

    double get_time() {
        get_time_calls++;
        timeval tv;
        gettimeofday(&tv, 0);
        return tv.tv_sec+tv.tv_usec*1e-6;
    }

    void print_time(string s) {
#ifdef LOCAL
    double x = get_time();
    fprintf(stderr,"%s cur=%.6lf lap=%.6lf\n",s.c_str(),x,x-last_call);
    last_call = x;
#endif
    }

    void init_time() {
        start_time = get_time();
        last_call = start_time;
    }
}
/******************************************************************/



/*************PARAMETER******************/
const int KMER = 32, PREFIX = 14, SHIFT = 2*PREFIX, SHIFTLEFT = 64 - (2*PREFIX);
uint64_t RIGHT31 = 4611686018427387903ULL;
uint32_t RIGHT16 = 4294967295, MIDDLE[8], MAXBIT = 1<<(2*PREFIX);
int nameFamily[3000005], nameGenus[3000005];
LL globalCom = 0;
/****************************************/

struct OutputData {
    uint64_t block_id;
    string dataString;
};


class HashTable {
private:
    uint32_t *stID, *suffix, *taxoID;
    size_t cntHash;

public:
    HashTable(){};
    ~HashTable(){
        delete[] suffix;
        delete[] taxoID;
        delete[] stID;
    };

    int distStringDP(uint32_t u, uint32_t v) {
        int dp[11][11];
        memset(dp, 0, sizeof(dp));
        int m = 10, n = 10;
        for (int i = 0; i <= m; i++)
            dp[i][0] = i;
        for (int j = 0; j <= n; j++)
            dp[0][j] = j;

        for (int i = 1; i <= m; i++) {
            int down = max(i-2,1), up = min(i+2, n);
            for (int j = down; j <= up; j++) {
                if(getbit(u,(i-1)*2) != getbit(v,(j-1)*2) || getbit(u,(i-1)*2+1) != getbit(v,(j-1)*2+1)){
                    dp[i][j] = min( 1 + dp[i-1][j],
                                   min(1 + dp[i][j-1],
                                    1 + dp[i-1][j-1]));
                }
                else
                    dp[i][j] = dp[i-1][j-1];
            }
        }
        return dp[m][n];
    }

    void init(uint64_t sz) {
        cntHash = 0;
        suffix = new uint32_t[sz+1];
        taxoID = new uint32_t[sz+1];
        stID = new uint32_t[MAXBIT+1];
    }

    II check_approximate(uint32_t id, uint32_t val) {
        int cnt = 0;
        FO (i, stID[id], stID[id+1]) {
            if ((suffix[i] & MIDDLE[1]) == (val & MIDDLE[1]) ||
                (suffix[i] & MIDDLE[2]) == (val & MIDDLE[2]) ||
                (suffix[i] & MIDDLE[3]) == (val & MIDDLE[3]) ||
                (suffix[i] & MIDDLE[4]) == (val & MIDDLE[4]) ||
                (suffix[i] & MIDDLE[5]) == (val & MIDDLE[5]) ||
                (suffix[i] & MIDDLE[6]) == (val & MIDDLE[6]) ) {
                    cnt++;
                    if (distStringDP(suffix[i], val) <= 2)
                        return II(taxoID[i], cnt);
               }
        }
        return II(0, cnt);
    }

    void read(string file) {
        string fileSize = file + "_Size";
        string fileSuffix = file + "_Suffix";
        string fileTaxo = file + "_Taxo";
        ifstream ifsSize(fileSize.c_str());
        ifstream ifsSuffix(fileSuffix.c_str());
        ifstream ifsTaxo(fileTaxo.c_str());

        uint64_t cntDB = 0;
        ifsSize.read((char *) &cntDB, sizeof(cntDB));
        DEBUG(cntDB);
        init(cntDB+1);

        FO (id,0,MAXBIT) {
            uint32_t num;
            ifsSize.read((char *) &num, sizeof(num));
            stID[id] = cntHash+1;
            uint32_t val[num], taxa[num];
            memset(val, 0, sizeof(val));
            memset(taxa, 0, sizeof(taxa));
            ifsSuffix.read((char *) &val, sizeof(val));
            ifsTaxo.read((char *) &taxa, sizeof(taxa));
            FO (i,0,num) {
                ++cntHash;
                suffix[cntHash] = val[i];
                taxoID[cntHash] = taxa[i];
            }
        }
        ifsSize.close(); ifsSuffix.close(); ifsTaxo.close();
        stID[MAXBIT] = cntHash;
    }

};
// Variable
HashTable HT;

inline int get_code(char c) {
    if (c == 'A') return 0;
    if (c == 'C') return 1;
    if (c == 'G') return 2;
    if (c == 'T') return 3;
    return 0;
}

uint64_t toNumDNA(string &s, int a, int len) {
    uint64_t ans = 0;
    for (int i = a; i < a+len; i++) {
        ans <<= 2;
        ans |= get_code(s[i]);
    }
    return ans;
}

inline uint64_t reverseMask(uint64_t _ikmer, int m_k) {
    uint64_t _ikmerR = _ikmer;
    // The following 6 lines come from Jellyfish source code
    _ikmerR = ((_ikmerR >> 2)  & 0x3333333333333333UL) | ((_ikmerR & 0x3333333333333333UL) << 2);
    _ikmerR = ((_ikmerR >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_ikmerR & 0x0F0F0F0F0F0F0F0FUL) << 4);
    _ikmerR = ((_ikmerR >> 8)  & 0x00FF00FF00FF00FFUL) | ((_ikmerR & 0x00FF00FF00FF00FFUL) << 8);
    _ikmerR = ((_ikmerR >> 16) & 0x0000FFFF0000FFFFUL) | ((_ikmerR & 0x0000FFFF0000FFFFUL) << 16);
    _ikmerR = ( _ikmerR >> 32                        ) | (_ikmerR                        << 32);
    _ikmerR = (((uint64_t)-1) - _ikmerR) >> (64 - (m_k << 1));
    return _ikmerR;
}

II ClassifySequence(string &s, HashTable &HT) {
    int lenSeq = s.size();
    if(lenSeq < 100)
        return II(-1, 0);

    vector<uint32_t> ans, Vid, Vval;
    uint64_t tt = toNumDNA(s, 0, KMER);
    uint64_t tmp = reverseMask(tt, KMER);
    if(tmp > tt)
        tmp = tt;
    Vid.push_back(tmp >> SHIFTLEFT);
    Vval.push_back(tmp & RIGHT16);

    FO (i,1,lenSeq-KMER) {
        tt = ((tt & RIGHT31) << 2) | get_code(s[i+31]);
        tmp = reverseMask(tt, KMER);
        if (tmp > tt)
            tmp = tt;
        Vid.push_back(tmp >> SHIFTLEFT);
        Vval.push_back(tmp & RIGHT16);
    }

    int readCom = 0;
    FO (i,0,Vid.size()) {
        II result_match = HT.check_approximate(Vid[i], Vval[i]);
        readCom += result_match.second;
        if (result_match.first > 0)
            ans.push_back(result_match.first);
    }
    if (ans.size() == 0)
        return II(-1, readCom);

    vector<uint32_t> VGenus, VSpecies;
    for (auto i : ans) {
        if (nameGenus[i] != 0)
            VGenus.push_back(nameGenus[i]);
        else
            VGenus.push_back(i);
    }

    sort(VGenus.begin(), VGenus.end());
    VGenus.push_back(0);
    int finalGenus = 0, freqGenus = 1, cntGenus = 0, maxx = 0;
    FO (i,0,VGenus.size()-1) {
        if (VGenus[i] == VGenus[i+1])
            freqGenus++;
        else {
            cntGenus++;
            if(freqGenus > maxx){
                maxx = freqGenus;
                freqGenus = 1;
                finalGenus = VGenus[i];
            }
        }
    }

    int cntTaxa = 1, finalTaxa = 0, cntHit = maxx;
    if (cntHit == 1 && (cntGenus >= 2 || lenSeq >= 500))
        finalTaxa = -1;
    else if(cntHit == 2 && cntGenus >= 4)
        finalTaxa = -1;
    else{
        for (auto i : ans) {
            if (nameGenus[i] == finalGenus)
                VSpecies.push_back(i);
        }
        sort(VSpecies.begin(), VSpecies.end());
        VSpecies.push_back(0);
        maxx = 0;
        finalTaxa = finalGenus;
        FO (i,0,VSpecies.size()-1) {
            if (VSpecies[i] == VSpecies[i+1])
                cntTaxa++;
            else{
                if(cntTaxa > maxx) {
                    maxx = cntTaxa;
                    cntTaxa = 1;
                    finalTaxa = VSpecies[i];
                }
            }
        }
    }
    return II(finalTaxa, readCom);
}

void usage(){
    cerr << "./CDKAM.sh DBname input output --fasta\n Or \n";
    cerr << "./CDKAM.sh DBname input output --fasta nthread N\n";
}

int main (int argc, char **argv) {
    if (argc != 6 && argc != 8) {
        usage();
        exit(1);
    }
    int num_threads = 1;
    if (argc == 8 && string(argv[6]) == "nthread")
        num_threads = atoi(argv[7]);

    omp_set_num_threads(num_threads);
    DEBUG(num_threads);


    ifstream finFGS(argv[2]);
    int valFamily, valGenus, valSpecies, valOrder;
    while (finFGS >> valFamily >> valGenus >> valSpecies ) {
        nameFamily[valSpecies] = valFamily;
        nameGenus[valSpecies] = valGenus;
    }
    finFGS.close();

    FOR (k,1,6) {
        MIDDLE[k] = 0;
        FO (i,10,16) if(i != 9+k) {
            MIDDLE[k] = onbit(MIDDLE[k], 2*i);
            MIDDLE[k] = onbit(MIDDLE[k], 2*i+1);
        }
    }


    // Read HashTable file
    string file(argv[1]);
    cerr << "Loading database ..." << endl;
    double main_time = Time::get_time();
    HT.read(file);
    double read_time = Time::get_time() - main_time;
    DEBUG(read_time);
    printRam();


    cerr << "Start testing" << endl;
    ifstream finReads(argv[3]);
    ofstream fout(argv[4]);
    string mode(argv[5]);

    int formatID = 0;
    if (mode == "--fasta") {
        DEBUG("FASTA");
        formatID = 0;
    }
    else if (mode == "--fastq") {
        DEBUG("FASTQ");
        formatID = 1;
    }

    // The priority queue for output is designed to ensure fragment data
    // is output in the same order it was input
    auto comparator = [](const OutputData &a, const OutputData &b) {
        return a.block_id > b.block_id;
    };
    std::priority_queue<OutputData, vector<OutputData>, decltype(comparator)>
    output_queue(comparator);
    uint64_t next_input_block_id = 0;
    uint64_t next_output_block_id = 0;
    omp_lock_t output_lock;
    omp_init_lock(&output_lock);

    #pragma omp parallel
    {
        Reader RR;
        RR.file_format_ = formatID;
        string headID, seq;
        uint64_t block_id;
        OutputData out_data;

        while(true) {
            int threadCom = 0;
            bool ok_read = false;
            #pragma omp critical(seqread)
            {
                ok_read = RR.LoadBlock(finReads, (size_t)3*1024*1024);
                block_id = next_input_block_id++;
            }

            if (! ok_read)
                break;

            // Reset ostringstream
            ostringstream oss;
            oss.str("");
            while (true) {
                auto valid_fragment = RR.NextSequence(headID, seq);
                if (! valid_fragment)
                    break;

                II ans = ClassifySequence(seq, HT);
                threadCom += ans.second;
                oss << headID << "\t" << seq.size() << "\t" << ans.first << "\n";
            }

            #pragma omp atomic
            globalCom += threadCom;

            out_data.block_id = block_id;
            out_data.dataString.assign(oss.str());

            #pragma omp critical(output_queue)
            {
                output_queue.push(out_data);
            }

            bool output_loop = true;
            while (output_loop) {
                #pragma omp critical(output_queue)
                {
                    output_loop = ! output_queue.empty();
                    if (output_loop) {
                        out_data = output_queue.top();
                        if (out_data.block_id == next_output_block_id) {
                            output_queue.pop();
                            // Acquiring output lock obligates thread to print out
                            // next output data block, contained in out_data
                            omp_set_lock(&output_lock);
                            next_output_block_id++;
                        }
                    else
                        output_loop = false;
                    }
                }
            if (! output_loop)
                break;
            // Past this point in loop, we know lock is set

            fout << out_data.dataString;
            omp_unset_lock(&output_lock);
            }
        }
    } // end parallel block

    omp_destroy_lock(&output_lock);

    double test_time = Time::get_time() - main_time;
    DEBUG(test_time);
    finReads.close();
    fout.close();


    return 0;
}


