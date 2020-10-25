/**
# CDKAM: a metagenomic classification tool using discriminative k-mers and approximate matching strategy
# Copyright 2019-2020
# Department of Bioinformatics and Biostatistics, Shanghai Jiao Tong University
# Contact information: buikien.dp@sjtu.edu.cn, ccwei@sjtu.edu.cn
#
# Function: Collecting kmers at species and genus rank level
*/

#include <bits/stdc++.h>
#include "helpers.h"
#define FOR(i,a,b) for(size_t i=a;i<=b;i++)
#define FO(i,a,b) for(size_t i=a;i<b;i++)
using namespace std;

#define PRINTLOG 1

const int KMER = 32, TWO22 = 4194304,  BIT11 = 4194303, RANGE = 7, LIMITgenus = 3;
const uint64_t LEFT31 = 4611686018427387903ULL;
/// RANGE = 5, X = 20%
/// RANGE = 7, X = 15%
/// RANGE = 10, X = 10%
/// RANGE = 20, X = 5%


class HashTable {
private:
    vector<uint64_t> mTable[TWO22+1], Vtmp;

public:
    HashTable(){};

    void init() {
        for (size_t i = 0; i < TWO22; ++i) {
            mTable[i].clear();
            mTable[i].shrink_to_fit();
        }
    }

    void insert(uint64_t num) {
        uint64_t id = num & BIT11;
        uint64_t val = num >> 22;
        mTable[id].push_back(val);
    }

    void sortData() {
        for (size_t i = 0; i < TWO22; ++i) {
            if(mTable[i].size() >= 3) {
                sort(mTable[i].begin(), mTable[i].end());
                Vtmp.clear();
                Vtmp.shrink_to_fit();
                Vtmp.push_back(mTable[i][0]);
                Vtmp.push_back(mTable[i][1]);
                for (size_t j = 2; j < mTable[i].size(); ++j) {
                    if (mTable[i][j] != mTable[i][j-2]) {
                        Vtmp.push_back(mTable[i][j]);
                    }
                }
                mTable[i].clear();
                mTable[i].shrink_to_fit();
                mTable[i] = Vtmp;
            }
        }
    }

    int check_freq(uint64_t num) {
        uint64_t id = num & BIT11;
        uint64_t val = num >> 22;
        size_t pos1 = lower_bound(mTable[id].begin(), mTable[id].end(), val) - mTable[id].begin();
        size_t pos2 = upper_bound(mTable[id].begin(), mTable[id].end(), val) - mTable[id].begin();
        // == 1 found it
        return pos2 - pos1;
    }
};


LL sumSet = 0, sumFinalKMER = 0, sumKMERGenome = 0;
set<uint64_t> S[TWO22+1];
vector<uint64_t> VGkmer, VPFamily;
vector<uint32_t> Vtaxa;
HashTable HT_Family, HT_Genus;


void convertUpperCase(string &s) {
    /// all are uppercase
    for (size_t i = 0; i < s.size(); ++i) {
        if(s[i] >= 'a' && s[i] <= 'z' && s[i] != 'x')
            s[i] += 'A'-'a';
    }
}

void insertKmerIntoTWO22Sets(uint64_t &tmp, const uint32_t speciesTID, const int flag) {
    if(tmp == 0) return;

    uint64_t tmp2 = reverseMask3(tmp, KMER);
    if(tmp > tmp2) tmp = tmp2;

    uint64_t setID, val;
    setID = tmp & BIT11;
    val = tmp >> 22;
    if (S[setID].find(val) == S[setID].end()) {
        S[setID].insert(val);
        VGkmer.push_back(tmp);
        Vtaxa.push_back(speciesTID);
        if(flag) VPFamily.push_back(tmp);
    }
}

void solveTaxa(const uint32_t speciesTID, const int stepID, vector<IIS> &vectorSpeciesList, const int flag) {
    LL numKmerOfSpecies = 0;
    for (size_t i = 0; i < TWO22; ++i) {
        if(!S[i].empty()) S[i].clear();
    }
    for (auto i : vectorSpeciesList) {
        string file = i.first;
        ifstream fin2(file.c_str());
        string genome = "", s, genID;

        while (getline(fin2, s)) {
            if (s[0] == '>' && genome.size() > 0) {
                convertUpperCase(genome);

                uint64_t tmp;
                for (size_t i = 0; i < genome.size() - KMER; ++i) {
                    tmp = toNumDNA(genome, i, KMER);
                    insertKmerIntoTWO22Sets(tmp, speciesTID, flag);
                }
                numKmerOfSpecies += genome.size() - KMER;

                genome = ""; /// reset size = 0
            }
            else genome += s;
        }

        {
                /// last sequence
                convertUpperCase(genome);

                uint64_t tmp;
                for (size_t i = 0; i < genome.size() - KMER; ++i) {
                    tmp = toNumDNA(genome, i, KMER);
                    insertKmerIntoTWO22Sets(tmp, speciesTID, flag);
                }
                numKmerOfSpecies += genome.size() - KMER;
        }
        fin2.close();
    }
    sumKMERGenome += numKmerOfSpecies;
    if (PRINTLOG)
        cerr << "step = " << stepID  << " SpeciesID = " << speciesTID << "  :  " << numKmerOfSpecies << "  " << VGkmer.size() << endl;
}

void solveGenus(const int familyTID, const int genusTID, const int stepID, vector<IIS> &vectorGenusList, ofstream &foutF, ofstream &ofs) {
    bool isFamilyFlag = false;
    if (stepID <= LIMITgenus) isFamilyFlag = true;

    map<int, int> mapTaxa;
    int cntTaxa = 0, iVtaxa[20000];
    map<int, int> freqSpecies;
    vector<II> filesInSpecies;
    vector<IIS> vectorSpecies[20000];
    VGkmer.clear(); VGkmer.shrink_to_fit();
    Vtaxa.clear(); Vtaxa.shrink_to_fit();
    HT_Genus.init();

    for (auto i : vectorGenusList) {
        int valSpecies = i.second;
        freqSpecies[valSpecies]++;
    }
    /// sorting vectorSpecies according to the size
    for (auto i : freqSpecies) {
        filesInSpecies.push_back(II(i.second, i.first));
    }
    sort(filesInSpecies.begin(), filesInSpecies.end(), greater<II> ());
    for (auto i : filesInSpecies) {
        int valSpecies = i.second;
        mapTaxa[valSpecies] = ++cntTaxa;
        iVtaxa[cntTaxa] = valSpecies;
    }

    for (auto i : vectorGenusList) {
        string file = i.first;
        int val = i.second;
        vectorSpecies[mapTaxa[val]].push_back(IIS(file, val));
    }


    /** FoutF print */
    FOR (i,1,cntTaxa) {
        foutF << familyTID << " " << genusTID  << " " << iVtaxa[i] << endl;
    }
    /** END */


    FOR (i,1,cntTaxa) {
        if (isFamilyFlag) solveTaxa(iVtaxa[i], i, vectorSpecies[i], 1);
        else              solveTaxa(iVtaxa[i], i, vectorSpecies[i], 0);
    }

    /**  HT_Genus */
    for (auto i : VGkmer) {
        HT_Genus.insert(i);
    }
    HT_Genus.sortData();
    /** END */


    int cntans = 0;
    if (isFamilyFlag) {
        for (size_t i = 0; i < VGkmer.size(); i += RANGE) {
            size_t pos = i;
            int ok = 0;
            for (size_t k = 0; k < RANGE; k++) {
                pos = i+k;
                if (pos >= VGkmer.size()) break;
                if (HT_Genus.check_freq(VGkmer[pos]) == 1) { ok = 1;  break; }
            }
            if (ok == 0) { pos = i; Vtaxa[pos] = genusTID; }

            ofs.write((char *) &VGkmer[pos], sizeof(VGkmer[pos]));
            ofs.write((char *) &Vtaxa[pos], sizeof(Vtaxa[pos]));
            cntans++;
        }
    }
    else {
        for (size_t i = 0; i < VGkmer.size(); i += RANGE) {
            size_t pos = i;
            int ok = 0;
            for (size_t k = 0; k < RANGE; k++) {
                pos = i+k;
                if (pos >= VGkmer.size()) break;
                if (HT_Genus.check_freq(VGkmer[pos]) == 1 && HT_Family.check_freq(VGkmer[pos]) == 0) { ok = 1; break; }
            }
            if (ok == 1) {
                ofs.write((char *) &VGkmer[pos], sizeof(VGkmer[pos]));
                ofs.write((char *) &Vtaxa[pos], sizeof(Vtaxa[pos]));
                cntans++;
            }
        }
    }
    sumFinalKMER += cntans;
    sumSet += VGkmer.size();


    if(PRINTLOG)
        cerr << "GENUS = " <<  genusTID << " have " << cntTaxa << " taxas.    Num of KMERs:  " << cntans << " . VPFamily size = " << VPFamily.size() << "\n\n";
    if(PRINTLOG) printRam();
}

void solveGenusHuman(uint32_t speciesTID, vector<string> &humanFileList, ofstream &foutF, ofstream &ofs) {
    LL numKmerOfSpecies = 0;
    VPFamily.clear(); VPFamily.shrink_to_fit();
    Vtaxa.clear(); Vtaxa.shrink_to_fit();
    DEBUG("HUMAN");
    foutF << "9604 9605 9606\n";

    int step = 0;
    for (auto file : humanFileList) {
        step++;
        for (size_t i = 0; i < TWO22; ++i) {
            if(!S[i].empty()) S[i].clear();
        }
        ifstream fin2(file.c_str());
        string genome = "", s, genID;

        while (getline(fin2, s)) {
            if (s[0] == '>' && genome.size() > 0) {
                convertUpperCase(genome);

                uint64_t tmp;
                for (size_t i = 0; i < genome.size() - KMER; ++i) {
                    tmp = toNumDNA(genome, i, KMER);
                    insertKmerIntoTWO22Sets(tmp, speciesTID, 0);
                }
                numKmerOfSpecies += genome.size() - KMER;

                /// reset size = 0
                genome = "";
            }
            else genome += s;
        }

        {
                /// last sequence
                convertUpperCase(genome);

                uint64_t tmp;
                for (size_t i = 0; i < genome.size() - KMER; ++i) {
                    tmp = toNumDNA(genome, i, KMER);
                    insertKmerIntoTWO22Sets(tmp, speciesTID, 0);
                }
                numKmerOfSpecies += genome.size() - KMER;
        }
        fin2.close();

        if (PRINTLOG)
            cerr << "Human id = " << step << "  :  " << numKmerOfSpecies << "  " << VPFamily.size() << endl;
    }


    int cntans = 0;
    for (size_t pos = 0; pos < VPFamily.size(); pos += RANGE) {
        ofs.write((char *) &VPFamily[pos], sizeof(VPFamily[pos]));
        ofs.write((char *) &Vtaxa[pos], sizeof(Vtaxa[pos]));
        cntans++;
    }

    printRam();
    sumKMERGenome += numKmerOfSpecies;
    sumFinalKMER += cntans;
    sumSet += VPFamily.size();
    if(PRINTLOG)
        cerr << "Human " <<  speciesTID << " have Num of KMERs:  " << cntans << " . VPFamily size = " << VPFamily.size() << "\n\n";
}

void solveFamily(int familyTID, vector<IISS> &vectorFamilyList, ofstream &foutF, ofstream &foutD) {
    map<int, int> mapGenus;
    map<int, int> freqGenus;
    vector<II> filesInGenus;
    vector<IIS> vectorGenus[20000];
    int cntGenus = 0, idGenus[20000];
    VPFamily.clear(); VPFamily.shrink_to_fit();
    HT_Family.init();

    if (familyTID == 9604) {
        cerr << "Human genome" << endl;
        int valSpecies = vectorFamilyList[0].second.second;
        vector<string> humanFileList;
        for(auto i : vectorFamilyList) {
            humanFileList.push_back(i.first);
        }
        solveGenusHuman(valSpecies, humanFileList, foutF, foutD);
        return;
    }

    for (auto i : vectorFamilyList) {
        int valGenus = i.second.first;
        int valSpecies = i.second.second;
        freqGenus[valGenus]++;
    }
    /// sorting vectorGenus according to the size
    for (auto i : freqGenus) {
        filesInGenus.push_back(II(i.second, i.first));
    }
    sort(filesInGenus.begin(), filesInGenus.end(), greater<II> ());
    for (auto i : filesInGenus) {
        int valGenus = i.second;
        mapGenus[valGenus] = ++cntGenus;
        idGenus[cntGenus] = valGenus;
    }

    for (auto i : vectorFamilyList) {
        string file = i.first;
        int valGenus = i.second.first;
        int valSpecies = i.second.second;
        vectorGenus[mapGenus[valGenus]].push_back(IIS(file, valSpecies));
    }

    int limit = min(LIMITgenus, cntGenus);
    FOR (i,1,limit) {
        if(PRINTLOG)
            cerr << filesInGenus[i-1].second << " GENUS have  " << filesInGenus[i-1].first  << endl;
        solveGenus(familyTID, idGenus[i], i, vectorGenus[i], foutF, foutD);
    }


    /// Insert LIMITgenus kmer into HT_Family
    for (auto i : VPFamily) {
        HT_Family.insert(i);
    }
    HT_Family.sortData();

    FOR (i,limit+1,cntGenus) {
        if(PRINTLOG)
            cerr << filesInGenus[i-1].second << " GENUS have  " << filesInGenus[i-1].first << endl;
        solveGenus(familyTID, idGenus[i], i, vectorGenus[i], foutF, foutD);
    }

    if(PRINTLOG) {
        cerr << "FAMILY = " << familyTID << " have " << vectorFamilyList.size() << endl;
        cerr << "************************************************************\n\n";
    }
}

void usage() {
    cerr << "./buildDB targets.txt nameFamily.txt database\n";
}

int main(int argc, char **argv) {
    if(argc != 4) { usage(); exit(1); }

    ifstream fin(argv[1]);
    ofstream foutF(argv[2]);
    ofstream foutD(argv[3], ofstream::binary);

    vector<IISS> vectorFamily[20000];
    map<int, int> mapFamily;
    int *idFamily, *parent;
    idFamily = new int [20000];
    parent = new int [3000005];
    memset(parent, 0, sizeof(parent));
    int maxx = 2700000;
    string file;
    int t[8], step = 0, cntFamily = 0;


    while (fin >> file) {
        step++;
        /// NEW: 7 levels
        FOR (i, 1, 7) fin >> t[i];
        if(t[1] != -1 ) {
            int orderID, familyID, genusID, speciesID;

            speciesID = t[2];
            if(t[3] != -1) { genusID = t[3]; parent[speciesID] = genusID; }
            else{
                if(parent[speciesID] == 0) genusID = ++maxx;
                else genusID = parent[speciesID];
            }

            if(t[4] != -1) { familyID = t[4]; parent[genusID] = familyID; }
            else{
                if(parent[genusID] == 0) familyID = ++maxx;
                else familyID = parent[genusID];
            }

            if(mapFamily[familyID] == 0) {
                mapFamily[familyID] = ++cntFamily;
                idFamily[cntFamily] = familyID;
            }
            vectorFamily[mapFamily[familyID]].push_back(IISS(file, II(genusID, speciesID )));
        }
    }

    if (RANGE == 5) cerr << "Program selects X = 20%" << endl;
    else if (RANGE == 7) cerr << "Program selects X = 15%" << endl;
    else if (RANGE == 10) cerr << "Program selects X = 10%" << endl;
    else if (RANGE == 20) cerr << "Program selects X = 5%" << endl;
    cerr << endl << maxx << endl;

    FOR (i, 1, cntFamily) {
        solveFamily(idFamily[i], vectorFamily[i], foutF, foutD);
    }

    DEBUG(sumKMERGenome);
    DEBUG(sumSet);
    DEBUG(sumFinalKMER);
    cerr << endl << endl;
    foutF.close();
    foutD.close();

    return 0;
}
