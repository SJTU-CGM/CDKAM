/**
# CDKAM: a metagenomic classification tool using discriminative k-mers and approximate matching strategy
# Copyright 2019-2020
# Department of Bioinformatics and Biostatistics, Shanghai Jiao Tong University
# Contact information: buikien.dp@sjtu.edu.cn, ccwei@sjtu.edu.cn
#
# Function: Removing kmers that are shared by at least two genomes
*/

#include <bits/stdc++.h>
#include "helpers.h"

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
        for (size_t i = 0; i < MAXBIT; ++i) {
            if (mTable[i].size() > 2) {
                sort(mTable[i].begin(), mTable[i].end());
                Vtmp.clear();
                Vtmp.shrink_to_fit();
                int freq = 1;
                mTable[i].push_back(II32(-1,-1));
                for (size_t j = 0; j < mTable[i].size()-1; ++j) {
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
    for (size_t i = 0; i < MAXBIT; ++i) {
        HTsize += HT.mTable[i].size();
    }
    DEBUG(HTsize);
    foutSize.write((char *) &HTsize, sizeof(HTsize));

	for (size_t i = 0; i < MAXBIT; ++i) {
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
