/**
# CDKAM: a metagenomic classification tool using discriminative k-mers and approximate matching strategy
# Copyright 2019-2020
# Department of Bioinformatics and Biostatistics, Shanghai Jiao Tong University
# Contact information: buikien.dp@sjtu.edu.cn, ccwei@sjtu.edu.cn
#
# Function: Removing kmers that are shared by at least two genomes
*/

#include <bits/stdc++.h>
#define LL long long
#define ULL unsigned long long
#define FOR(i,a,b) for(uint32_t i=a;i<=b;i++)
#define FO(i,a,b) for(uint32_t i=a;i<b;i++)
#define DEBUG(a) {cerr << #a << ": " << (a) << endl; fflush(stderr); }

using namespace std;
template<class T> int getbit(T s, int i) { return (s >> i) & 1; }
template<class T> T onbit(T s, int i) { return s | (T(1) << i); }
template<class T> T offbit(T s, int i) { return s & (~(T(1) << i)); }
template<class T> int cntbit(T s) { return __builtin_popcount(s);}
typedef pair<uint32_t, uint32_t> II32;


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




inline int get_code(char c) {
    if (c == 'A') return 0;
    if (c == 'C') return 1;
    if (c == 'G') return 2;
    if (c == 'T') return 3;
    return 0;
}

inline char print_code(int c) {
    if (c == 0) return 'A';
    if (c == 1) return 'C';
    if (c == 2) return 'G';
    if (c == 3) return 'T';
    return 'N';
}





const int KMER = 32, PREFIX = 14, MAXBIT = 1<<(2*PREFIX);
uint64_t RIGHT16 = 4294967295;


class HashTable {
public:
    vector<II32> *mTable, Vtmp;

    HashTable(){};
    ~HashTable(){
        delete[] mTable;
    };

    void init(){
        mTable = new vector<II32> [MAXBIT];
    }

    void insert(int id, uint32_t val, uint32_t taxa) {
        mTable[id].push_back(II32(val, taxa));
    }

    void sortData() {
        int cntRemove = 0;
        FO (i,0,MAXBIT) {
            if (mTable[i].size() > 2) {
                sort(mTable[i].begin(), mTable[i].end());
                Vtmp.clear();
                Vtmp.shrink_to_fit();
                int freq = 1;
                mTable[i].push_back(II32(-1,-1));
                FO (j,0,mTable[i].size()-1) {
                    if (mTable[i][j].first != mTable[i][j+1].first) {
                        if(freq == 1)
                            Vtmp.push_back(mTable[i][j]);
                        else {
                            freq = 1;
                            cntRemove++;
                        }
                    }
                    else {
                        freq++;
                        cntRemove++;
                    }
                }
                mTable[i].clear();
                mTable[i].shrink_to_fit();
                mTable[i] = Vtmp;
            }
        }
        DEBUG(cntRemove);
    }

};
HashTable HT;


uint64_t reverseMask3(uint64_t _ikmer, int m_k) {
    uint64_t _ikmerR = _ikmer;
    // The following 6 lines come from Jellyfish source code
    _ikmerR = ((_ikmerR >> 2)  & 0x3333333333333333UL) | ((_ikmerR & 0x3333333333333333UL) << 2);
    _ikmerR = ((_ikmerR >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_ikmerR & 0x0F0F0F0F0F0F0F0FUL) << 4);
    _ikmerR = ((_ikmerR >> 8)  & 0x00FF00FF00FF00FFUL) | ((_ikmerR & 0x00FF00FF00FF00FFUL) << 8);
    _ikmerR = ((_ikmerR >> 16) & 0x0000FFFF0000FFFFUL) | ((_ikmerR & 0x0000FFFF0000FFFFUL) << 16);
    _ikmerR = ( _ikmerR >> 32                        ) | (_ikmerR                        << 32);
    _ikmerR = (((uint64_t)-1) - _ikmerR) >> (64 - (m_k << 1));
    if(_ikmer < _ikmerR) return _ikmer;
    return _ikmerR;
}

void usage() {
    cerr << "./DTB_unique database.txt FinalDTB\n";
}

int main(int argc, char **argv) {
    if(argc != 3) { usage(); exit(1); }

    ifstream ifs(argv[1]);
    string file(argv[2]);
    string fileSize = file + "_Size";
    string fileSuffix = file + "_Suffix";
    string fileTaxo = file + "_Taxo";
    ofstream foutSize(fileSize.c_str(), ofstream::binary);
    ofstream foutSuffix(fileSuffix.c_str(), ofstream::binary);
    ofstream foutTaxo(fileTaxo.c_str(), ofstream::binary);

    uint64_t cntDB = 0;
    HT.init();
    double main_time = Time::get_time();
    DEBUG("Start reading");

    uint32_t taxaID;
    uint64_t tmp;
    int SHIFTLEFT = 64 - (2*PREFIX);
    while (ifs.read((char *) &tmp, sizeof(tmp))) {
        ifs.read((char *) &taxaID, sizeof(taxaID));
        cntDB++;
        uint64_t id = tmp >> SHIFTLEFT; /// PREFIX bits left
        uint64_t val = tmp & RIGHT16;
        HT.insert(id, (uint32_t) val, taxaID);
	}
	DEBUG(cntDB);
	ifs.close();

	double read_time = Time::get_time() - main_time;
    DEBUG(read_time);

	HT.sortData();
    double sort_time = Time::get_time() - main_time;
    DEBUG(sort_time);


    DEBUG("Writing the final database");
    uint64_t cntNum = 0, HTsize = 0;
    FO (i,0,MAXBIT) {
        HTsize += HT.mTable[i].size();
    }
    DEBUG(HTsize);
    foutSize.write((char *) &HTsize, sizeof(HTsize));

	FO (i,0,MAXBIT) {
	    uint32_t sz = HT.mTable[i].size();
	    foutSize.write((char *) &sz, sizeof(sz));
	    if (sz == 0)
            cntNum++;
        for (auto u : HT.mTable[i]) {
            //cout << u.fi << " " << u.se << "\n";
            foutSuffix.write((char *) &u.first, sizeof(u.first));
            foutTaxo.write((char *) &u.second, sizeof(u.second));
        }
	}
    foutSize.close();
    foutSuffix.close();
    foutTaxo.close();

    return 0;
}
