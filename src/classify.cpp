/**
# CDKAM: a metagenomic classification tool using discriminative k-mers and approximate matching strategy
# Copyright 2019-2020
# Department of Bioinformatics and Biostatistics, Shanghai Jiao Tong University
# Contact information: buikien.dp@sjtu.edu.cn, ccwei@sjtu.edu.cn
#
# Function: Reading sequences and carrying out classification
*/



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
const int KMER = 32, PREFIX = 14, SHIFT = 2* PREFIX, SHIFTLEFT = 64 - (2*PREFIX), MAXBIT = 1<<(2*PREFIX);
uint64_t RIGHT31 = 4611686018427387903ULL;
uint32_t  RIGHT10 = 1048575, RIGHT16 = 4294967295, MIDDLE[8];
int nameFamily[3000005], nameGenus[3000005];
/***************************************************/



class HashTable {
private:
    int dp[20][20];
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
        int m = 10, n = 10;
        for (int i = 0; i <= m; i++)
            dp[i][0] = i;
        for (int j = 0; j <= n; j++)
            dp[0][j] = j;

        for (int i = 1; i <= m; i++) {
            int down = max(i-2,1), up = min(i+2, n);
            for (int j = down; j <= up; j++) {
                if(getbit(u,(i-1)*2) != getbit(v,(j-1)*2) || getbit(u,(i-1)*2+1) != getbit(v,(j-1)*2+1)){
                ///if (a[i-1] != b[j-1]) {
                    dp[i][j] = min( 1 + dp[i-1][j],  // deletion
                                   min(1 + dp[i][j-1],  // insertion
                                    1 + dp[i-1][j-1])); // replacement
                }
                else  dp[i][j] = dp[i-1][j-1];
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

    int check_approximate(uint32_t  id, uint32_t val) {
        FO (i, stID[id], stID[id+1]) {
            if ((suffix[i] & MIDDLE[1]) == (val & MIDDLE[1]) ||
               (suffix[i] & MIDDLE[2]) == (val & MIDDLE[2]) ||
               (suffix[i] & MIDDLE[3]) == (val & MIDDLE[3]) ||
               (suffix[i] & MIDDLE[4]) == (val & MIDDLE[4]) ||
               (suffix[i] & MIDDLE[5]) == (val & MIDDLE[5]) ||
               (suffix[i] & MIDDLE[6]) == (val & MIDDLE[6]) )
                    if (distStringDP(suffix[i], val) <= 2)
                        return taxoID[i];
        }
        return 0;
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
            //FO(i,0,num) ifsSuffix.read((char *) &val[i], sizeof(val[i]));
            //FO(i,0,num) ifsSuffix.read((char *) &taxa[i], sizeof(taxa[i]));
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
/// Variable
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

int ClassifySequence(int step, string &s, HashTable &HT) {
    int lenSeq = s.size();
    if(lenSeq < 100)
        return -1;

    vector<uint32_t> ans, Vid, Vval;
    uint64_t tt = toNumDNA(s, 0, KMER);
    uint64_t tmp = reverseMask(tt, KMER);
    if(tmp > tt) tmp = tt;
    Vid.push_back(tmp >> SHIFTLEFT);
    Vval.push_back(tmp & RIGHT16);
    FO (i,1,lenSeq-KMER)   {
        tt = ((tt & RIGHT31) << 2) | get_code(s[i+31]);
        tmp = reverseMask(tt, KMER);
        if(tmp > tt) tmp = tt;
        Vid.push_back(tmp >> SHIFTLEFT);
        Vval.push_back(tmp & RIGHT16);
    }

    FO (i,0,Vid.size()) {
        int ok = HT.check_approximate(Vid[i], Vval[i]);
        if (ok > 0)
            ans.push_back(ok);
    }
    if (ans.size() == 0)
        return -1;

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
    if (cntHit == 1 && (cntGenus >= 3 || lenSeq >= 800))
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
    return finalTaxa;
}

void usage(){
    cerr << "./CDKAM database nameFamily input output --fasta\n";
}

int main(int argc, char **argv) {
    if (argc != 6) { usage(); exit(1); }

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



    /// Read HashTable file
    string file(argv[1]);
    DEBUG("Reading");
    double main_time = Time::get_time();
    HT.read(file);
    double read_time = Time::get_time() - main_time;
    DEBUG(read_time);


    DEBUG("Testing");
    printRam();
    ifstream finReads(argv[3]);
    ofstream fout(argv[4]);
    string mode(argv[5]);
    Reader RR;
    if (mode == "--fasta") {
        DEBUG("FASTA");
        RR.file_format_ = 0;
    }
    else if (mode == "--fastq") {
        DEBUG("FASTQ");
        RR.file_format_ = 1;
    }

    {
        int step = 0;
        string seq;
        vector<II> V;
        while(true) {
            bool ok_read = RR.LoadBlock(finReads, (size_t)3*1024*1024);
            if (! ok_read)
                break;
            while (true) {
                auto valid_fragment = RR.NextSequence(seq);
                if (! valid_fragment)
                    break;

                step++;
                int ans = ClassifySequence(step, seq, HT);
                fout << step << "\t" << seq.size() << "\t" << ans << "\n";
            }
        }
    }
    double test_time = Time::get_time() - main_time;
    DEBUG(test_time);
    finReads.close();
    fout.close();
    //*/

    return 0;
}


