/**
# CDKAM: a metagenomic classification tool using discriminative k-mers and approximate matching strategy
# Copyright 2019-2020
# Department of Bioinformatics and Biostatistics, Shanghai Jiao Tong University
# Contact information: buikien.dp@sjtu.edu.cn, ccwei@sjtu.edu.cn
#
# Function: Generating each species's genome from .fna file
*/

#include <bits/stdc++.h>
#define LL long long
#define ULL unsigned long long
#define FOR(i,a,b) for(int i=a;i<=b;i++)
#define FO(i,a,b) for(int i=a;i<b;i++)
#define DEBUG(a) {cerr << #a << ": " << (a) << endl; fflush(stderr); }

using namespace std;

typedef pair<int, int> II;
typedef pair<II, int> III;
typedef pair<string, int> IIS;
typedef pair<string, II> IISS;
#define maxn 3000005



string mapSpecies[maxn], mapGenus[maxn], mapFamily[maxn], mapOrder[maxn], mapClass[maxn], mapPhylum[maxn];
IIS parent[3000001];



vector<string> tokenize(const string& row) {
    vector<string> tokens;
    for(int i=0, pos=0, n=row.size(); i<n; ++i) {
        if(i==n-1 || row[i+1]=='|' || row[i+1] == ' ') {
            string token = row.substr(pos, (i+1)-pos);
            tokens.push_back(token);
            pos = i+2;
        }
    }
    return tokens;
}

inline int get_code(char c) {
    if (c == 'A') return 0;
    if (c == 'T') return 1;
    if (c == 'C') return 2;
    if (c == 'G') return 3;
    return 0;
}

inline char print_code(int c) {
    if (c == 0) return 'A';
    if (c == 1) return 'T';
    if (c == 2) return 'C';
    if (c == 3) return 'G';
    return 'N';
}

int toNum(string s) {
    int ans = 0;
    FO (i,0,s.size()) ans = ans*10 + s[i] - 48;
    return ans;
}

string toString(int n) {
    if (n == -1) return "-1";
    string ans = "";
    while (n){
        ans = char(n%10 + 48) + ans;
        n /= 10;
    }
    return ans;
}

uint64_t toNumDNA(string &s, int a, int len) {
    LL ans = 0;
    for (int i = a; i < a+len; i++) {
        ans <<= 2;
        ans |= get_code(s[i]);
    }
    return ans;
}


void usage(){
    cerr << "./genSpecies premap.txt nodes.dmp library.fna target.txt\n";
}


int main (int argc, char **argv) {
    if (argc != 5) {
        usage();
        exit(1);
    }
    ifstream fin(argv[1]);
    ifstream fin2(argv[2]);
    ifstream fin3(argv[3]);

    /// Read Taxonomy NCBI file
    FOR(i,1,3000000) parent[i] = IIS("",0);
    string t;
    while (getline(fin2, t)) {
        vector<string> V = tokenize(t);
        FOR(i,0,2)
            V[i].erase(V[i].size()-1,1);
        FOR(i,1,2)
            V[i].erase(0, 1);
        int u = toNum(V[0]);
        int v = toNum(V[1]);
        parent[u] = IIS(V[2], v);
    }
    fin2.close();

    /// Read Premap.txt
    string taxid;
	while (fin >> taxid) {
	    vector<string> V = tokenize(taxid); /// CDKAM|GCF_000091045.1|214684|NC_006670.1
        int strainID = toNum(V[2]);
        int phylumID = -1, classID = -1, orderID = -1, familyID = -1, genusID = -1, speciesID = strainID;
        int st = strainID;
        while (st != 131567) {
            if (st <= 2) break;
            if (parent[st].first == "species")      speciesID = st;
            else if(parent[st].first == "genus")    genusID = st;
            else if (parent[st].first == "family")  familyID = st;
            else if (parent[st].first == "order")   orderID = st;
            else if (parent[st].first == "class")   classID = st;
            else if (parent[st].first == "phylum")  phylumID = st;
            st = parent[st].second;
        }

        mapSpecies[strainID] = toString(speciesID);
        mapGenus[strainID] = toString(genusID);
        mapFamily[strainID] = toString(familyID);
        mapOrder[strainID] = toString(orderID);
        mapClass[strainID] = toString(classID);
        mapPhylum[strainID] = toString(phylumID);
        //cout << genID << " " << speciesID << " " << genusID << " " << familyID << " " << orderID << " " << classID << " " << phylumID << endl;
	}
    fin.close();

    /// Read library.fna and then split it into TaxaID files
    string s, oldFile = "";
    vector<string> Genome;
    vector<IIS> Vlink;
    while (getline(fin3, s)) {
        if (s[0] == '>') {
            if (oldFile.size() > 0) {
                vector<string> V = tokenize(oldFile);
                string file = "references/" + V[1] + "|" + V[2] + "|" + V[3] + ".txt";
                int id = toNum(V[2]);
                string link = string(argv[4]);
                link += "/" + file;

                ofstream ofs(link.c_str(), std::ios_base::app);
                ofs << oldFile << endl;
                for(auto i : Genome)
                    ofs << i << "\n";
                ofs.close();
                Genome.clear();
                Genome.shrink_to_fit();
                Vlink.push_back(IIS(link, id));
            }
            oldFile = s;
        }
        else{
            Genome.push_back(s);
        }
    }

    /// last genome
    {
                vector<string> V = tokenize(oldFile);
                string file = "references/" + V[1] + "|" + V[2] + "|" + V[3] + ".txt";
                int id = toNum(V[2]);
                string link = string(argv[4]);
                link += "/" + file;

                ofstream ofs(link.c_str(), std::ios_base::app);
                ofs << oldFile << endl;
                for(auto i : Genome)
                    ofs << i << "\n";
                ofs.close();
                Genome.clear();
                Genome.shrink_to_fit();
                Vlink.push_back(IIS(link, id));
    }
    fin3.close();

    /// Print library_name.txt, consists of file link and the full taxonomy path
    string output(argv[4]);
    output += ".txt";
    ofstream fout(output.c_str());
    sort(Vlink.begin(), Vlink.end());
    Vlink.erase(unique(Vlink.begin(), Vlink.end()), Vlink.end());
    for (auto i : Vlink) {
        int id = i.second;
        i.first += "\t" + toString(id) + "\t" + mapSpecies[id] + "\t" + mapGenus[id] + "\t" + mapFamily[id] + "\t" + mapOrder[id] + "\t" + mapClass[id] + "\t" + mapPhylum[id];
        fout << i.first << "\n";
    }
    fout.close();


	return 0;
}
