/**
# CDKAM: a metagenomic classification tool using discriminative k-mers and approximate matching strategy
# Copyright 2019-2020
# Department of Bioinformatics and Biostatistics, Shanghai Jiao Tong University
# Contact information: buikien.dp@sjtu.edu.cn, ccwei@sjtu.edu.cn
#
# Function: Collecting kmers at species and genus rank level
*/

#include <bits/stdc++.h>
#define LL long long
#define ULL unsigned long long
#define FOR(i,a,b) for(size_t i=a;i<=b;i++)
#define FO(i,a,b) for(size_t i=a;i<b;i++)
#define DEBUG(a) {cerr << #a << ": " << (a) << endl; fflush(stderr); }

using namespace std;

typedef pair<int, int> II;
typedef pair<string, int> IIS;
typedef pair<string, II> IISS;
#define PRINTLOG 0

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
        FO (i, 0, TWO22) {
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
        FO (i, 0, TWO22) if(mTable[i].size() >= 3) {
            sort(mTable[i].begin(), mTable[i].end());
            Vtmp.clear();
            Vtmp.shrink_to_fit();
            Vtmp.push_back(mTable[i][0]);
            Vtmp.push_back(mTable[i][1]);
            for (size_t j = 2; j < mTable[i].size(); j++) {
                if (mTable[i][j] != mTable[i][j-2]) {
                    Vtmp.push_back(mTable[i][j]);
                }
            }
            mTable[i].clear();
            mTable[i].shrink_to_fit();
            mTable[i] = Vtmp;
        }
    }

    int check_freq(uint64_t num) {
        uint64_t id = num & BIT11;
        uint64_t val = num >> 22;
        size_t pos1 = lower_bound(mTable[id].begin(), mTable[id].end(), val) - mTable[id].begin();
        size_t pos2 = upper_bound(mTable[id].begin(), mTable[id].end(), val) - mTable[id].begin();
        /// == 1 found it
        return pos2 - pos1;
    }
};


LL sumSet = 0, sumFinalKMER = 0, sumKMERGenome = 0;
set<uint64_t> S[TWO22+1];
vector<uint64_t> VGkmer, VPFamily;
vector<uint32_t> Vtaxa;
HashTable HT_Family, HT_Genus;


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

int toNum(string s) {
    int ans = 0;
    for(int i = 0; i < s.size(); i++)
        ans = ans*10 + s[i] - 48;
    return ans;
}

uint64_t toNumDNA(string &s, int a, int len) {
    uint64_t ans = 0;
    for (int i = a; i < a+len; i++) {
        if(s[i] == 'x') return 0;
        ans <<= 2;
        ans |= get_code(s[i]);
    }
    return ans;
}

inline uint64_t reverseMask3(uint64_t _ikmer, int m_k) {
    uint64_t _ikmerR = _ikmer;
    // m_k = size of Kmer
    // The following 6 lines come from Jellyfish source code
    _ikmerR = ((_ikmerR >> 2)  & 0x3333333333333333UL) | ((_ikmerR & 0x3333333333333333UL) << 2);
    _ikmerR = ((_ikmerR >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_ikmerR & 0x0F0F0F0F0F0F0F0FUL) << 4);
    _ikmerR = ((_ikmerR >> 8)  & 0x00FF00FF00FF00FFUL) | ((_ikmerR & 0x00FF00FF00FF00FFUL) << 8);
    _ikmerR = ((_ikmerR >> 16) & 0x0000FFFF0000FFFFUL) | ((_ikmerR & 0x0000FFFF0000FFFFUL) << 16);
    _ikmerR = ( _ikmerR >> 32                        ) | (_ikmerR                        << 32);
    _ikmerR = (((uint64_t)-1) - _ikmerR) >> (64 - (m_k << 1));
    return _ikmerR;
}


void solveTaxa(uint32_t speciesID, int id, vector<IIS> &vectorSpeciesList, int flag) {
    LL numKmerOfSpecies = 0;
    FOR (i, 0, TWO22) {
        if(!S[i].empty()) S[i].clear();
    }
    for (auto i : vectorSpeciesList) {
        string file = i.first;
        ifstream fin2(file.c_str());
        string genome = "", s, genID;

        while (getline(fin2, s)) {
            if (s[0] == '>' && genome.size() > 0) {
                /// all are uppercase
                FO (i,0,genome.size()) {
                    if(genome[i] >= 'a' && genome[i] <= 'z' && genome[i] != 'x')
                        genome[i] += 'A'-'a';
                }
                uint64_t group, tmp, id, val;
                FO (i, 0, genome.size() - KMER) {
                    tmp = toNumDNA(genome, i, KMER);
                    if (tmp == 0) continue;
                    uint64_t tmp2 = reverseMask3(tmp, KMER);
                    if(tmp > tmp2) tmp = tmp2;

                    id = tmp & BIT11;
                    val = tmp >> 22;
                    if (S[id].find(val) == S[id].end()) {
                        S[id].insert(val);
                        VGkmer.push_back(tmp);
                        Vtaxa.push_back(speciesID);
                        if(flag)
                            VPFamily.push_back(tmp);
                    }
                }
                numKmerOfSpecies += genome.size() - KMER;

                /// reset size = 0
                genome = "";
            }
            else genome += s;
        }

        {
                /// last sequence
                FO (i, 0, genome.size()) {
                    if(genome[i] >= 'a' && genome[i] <= 'z' && genome[i] != 'x')
                        genome[i] += 'A'-'a';
                }
                uint64_t group, tmp, id, val;
                FO (i, 0, genome.size() - KMER) {
                    tmp = toNumDNA(genome, i, KMER);
                    if(tmp == 0) continue;
                    uint64_t tmp2 = reverseMask3(tmp, KMER);
                    if(tmp > tmp2) tmp = tmp2;

                    id = tmp & BIT11;
                    val = tmp >> 22;
                    if (S[id].find(val) == S[id].end()) {
                        S[id].insert(val);
                        VGkmer.push_back(tmp);
                        Vtaxa.push_back(speciesID);
                        if(flag)
                            VPFamily.push_back(tmp);
                    }
                }
                numKmerOfSpecies += genome.size() - KMER;
        }
        fin2.close();
    }
    sumKMERGenome += numKmerOfSpecies;
    if (PRINTLOG)
        cerr << "step = " << id  << " SpeciesID = " << speciesID << "  :  " << numKmerOfSpecies << "  " << VGkmer.size() << endl;
}

void solveGenus(int Family, int Genus, int id, vector<IIS> &vectorGenusList, ofstream &foutF, ofstream &ofs) {
    map<int, int> mapTaxa;
    int cntTaxa = 0, iVtaxa[20000];
    map<int, int> freqSpecies;
    vector<II> Vtmp;
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
        Vtmp.push_back(II(i.second, i.first));
    }
    sort(Vtmp.begin(), Vtmp.end(), greater<II> ());
    for (auto i : Vtmp) {
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
        foutF << Family << " " << Genus  << " " << iVtaxa[i] << endl;
    }
    /** END */

    FOR (i,1,cntTaxa) {
        if (id <= LIMITgenus)
            solveTaxa(iVtaxa[i], i, vectorSpecies[i], 1);
        else
            solveTaxa(iVtaxa[i], i, vectorSpecies[i], 0);
    }

    /**  HT_Genus */
    for (auto i : VGkmer) {
        HT_Genus.insert(i);
    }
    HT_Genus.sortData();
    /** END */

    int cntans = 0;
    if (id <= LIMITgenus) {
        for (size_t i = 0; i < VGkmer.size(); i += RANGE) {
            size_t pos = i;
            int ok = 0;
            for (size_t k = 0; k < RANGE; k++) {
                pos = i+k;
                if (pos >= VGkmer.size()) break;
                if (HT_Genus.check_freq(VGkmer[pos]) == 1) { ok = 1;  break; }
            }
            if (ok == 0) { pos = i; Vtaxa[pos] = Genus; }

            ofs.write((char *) &VGkmer[pos], sizeof(VGkmer[pos]));
            ofs.write((char *) &Vtaxa[pos], sizeof(Vtaxa[pos]));
            cntans++;
        }
    }
    else{
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
        cerr << "GENUS = " <<  Genus << " have " << cntTaxa << " taxas.    Num of KMERs:  " << cntans << " . VPFamily size = " << VPFamily.size() << "\n\n";
}



void solveFamily(int Family, int id, vector<IISS> &vectorFamilyList, ofstream &foutF, ofstream &foutD) {
    map<int, int> mapGenus;
    map<int, int> freqGenus;
    vector<II> Vtmp;
    vector<IIS> vectorGenus[20000];
    int cntGenus = 0, idGenus[20000];
    VPFamily.clear(); VPFamily.shrink_to_fit();
    HT_Family.init();

    for (auto i : vectorFamilyList) {
        int valGenus = i.second.first;
        int valSpecies = i.second.second;
        freqGenus[valGenus]++;
    }
    /// sorting vectorGenus according to the size
    for (auto i : freqGenus) {
        Vtmp.push_back(II(i.second, i.first));
    }
    sort(Vtmp.begin(), Vtmp.end(), greater<II> ());
    for (auto i : Vtmp) {
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
            cerr << Vtmp[i-1].second << " GENUS have  " << Vtmp[i-1].first  << endl;
        solveGenus(Family, idGenus[i], i, vectorGenus[i], foutF, foutD);
    }

    for (auto i : VPFamily) {
        HT_Family.insert(i);
    }
    HT_Family.sortData();
    FOR (i,limit+1,cntGenus) {
        if(PRINTLOG)
            cerr << Vtmp[i-1].second << " GENUS have  " << Vtmp[i-1].first << endl;
        solveGenus(Family, idGenus[i], i, vectorGenus[i], foutF, foutD);
    }

    if(PRINTLOG) {
        cerr << "FAMILY = " << Family << " have " << vectorFamilyList.size() << endl;
        cerr << "************************************************************\n\n";
    }
}

void usage() {
    cerr << "./compressNEW targets.txt nameFamily.txt database\n";
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
    int t[8];
	int step = 0, cntFamily = 0;


	while (fin >> file) {
        step++;
		FOR (i, 1, 6) fin >> t[i];
		if(t[1] != -1 ) {
            int orderID, familyID, genusID, speciesID;

            speciesID = t[1];
            if(t[2] != -1) { genusID = t[2]; parent[speciesID] = genusID; }
            else{
                if(parent[speciesID] == 0) genusID = ++maxx;
                else genusID = parent[speciesID];
            }

            if(t[3] != -1) { familyID = t[3]; parent[genusID] = familyID; }
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

	DEBUG(maxx);
	DEBUG(RANGE);
	FOR (i, 1, cntFamily) {
        solveFamily(idFamily[i], i, vectorFamily[i], foutF, foutD);
	}
	DEBUG(sumKMERGenome);
    DEBUG(sumSet);
    DEBUG(sumFinalKMER);
    foutF.close();
    foutD.close();

	return 0;
}
