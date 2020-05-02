/**
# CDKAM: a metagenomic classification tool using discriminative k-mers and approximate matching strategy
# Copyright 2019-2020
# Department of Bioinformatics and Biostatistics, Shanghai Jiao Tong University
# Contact information: buikien.dp@sjtu.edu.cn, ccwei@sjtu.edu.cn
#
# Function: Translating result to Scientific Name and getting the Abundance
*/

#include <bits/stdc++.h>
#define LL long long
#define ULL unsigned long long
#define fi first
#define se second
#define FOR(i,a,b) for(int i=a;i<=b;i++)
#define FO(i,a,b) for(int i=a;i<b;i++)
#define DEBUG(a) {cerr << #a << ": " << (a) << endl; fflush(stderr); }

using namespace std;
typedef pair<LL, LL> II;
typedef pair<string, int> IIS;



vector<II> V;
vector<int> adj[3000005];
int testSpecies[3000005], testGenus[3000005], freq[3000005];
int nameFamily[3000005], nameGenus[3000005];
string MName[3000005];
IIS parent[3000005];


vector<string> tokenizeNode(const string& row) {
    vector<string> tokens;
    for(int i=0, pos=0, n=row.size(); i<n; ++i) {
        if(i==n-1 || row[i+1]=='|') {
            string token = row.substr(pos, (i+1)-pos);
            tokens.push_back(token);
            pos = i+2;
        }
    }
    return tokens;
}

int toNum(string s) {
    if(s == "NA" || s == "UNKNOWN") return -1;
    int ans = 0;
    FO(i,0,s.size()) ans = ans*10 + s[i] - 48;
    return ans;
}

int findGenus(int taxa) {
    if(taxa == -1) return -1;
    int st = taxa;
    while(st != 131567){
        if(st <= 2 ) break;
        if(parent[st].fi == "genus") return st;
        st = parent[st].se;
    }
    return taxa;
}

void usage(){
    cerr << "./translate DTB input output\n";
}


int main(int argc, char **argv) {
    if(argc != 4) { usage(); exit(1); }

    string DTB(argv[1]);
    string namesFile = DTB + "/taxonomy/names.dmp";
    string nodesFile = DTB + "/taxonomy/nodes.dmp";
    ifstream finNames(namesFile);
	ifstream finNodes(nodesFile);
	ifstream fin(argv[2]);
	ofstream fout(argv[3]);

    string s;
    while(getline(finNames, s)) {
        vector<string> V = tokenizeNode(s);
        V[0].erase(V[0].size()-1,1);
        V[1].erase(0,1);
        if(V[1].size() > 1)
            V[1].erase(V[1].size()-1,1);
        V[3].erase(0, 1);
        V[3].erase(V[3].size()-1,1);

        int ID = toNum(V[0]);
        if(MName[ID] == "" && V[3] == "scientific name"){
            MName[ID] = V[1];
        }
    }
    DEBUG("NAMES");

    FOR(i,1,3000000) parent[i] = IIS("",0);
    string t;
    while (getline(finNodes, t)) {
        vector<string> V = tokenizeNode(t);
        FOR (i,0,2) V[i].erase(V[i].size()-1,1);
        FOR (i,1,2) V[i].erase(0, 1);
        int u = toNum(V[0]);
        int v = toNum(V[1]);
        if(V[2] == "species") adj[v].push_back(u);
        parent[u] = IIS(V[2], v);
    }
    DEBUG("NODES");
    finNames.close();
    finNodes.close();

    DEBUG("TRANSLATION");
    string genID, id, seq;
    int stt, taxa, len;
	int step = 0, cntU = 0;
	while (fin >> stt >> len >> taxa) {
        testSpecies[++step] = taxa;
        testGenus[step] = findGenus(taxa);
        string ScientificName = "";
        if (testGenus[step] != -1)
            ScientificName = MName[testGenus[step]];
        fout << stt << "\t" << taxa << "\t" << testGenus[step];

        if (taxa == -1) {
            cntU++;
            fout << endl;
            continue;
        }
        if (taxa == testGenus[step]) {
            freq[taxa]++;
            fout << "\t (G) " << ScientificName << endl;
        }
        else {
            freq[taxa]++;
            freq[testGenus[step]]++;
            fout << "\t (S) " << ScientificName << endl;
        }
	}
	fout.close();


	ofstream foutA("abundance.txt");
	int cntReads = step, cntC = cntReads - cntU;
    vector<II> ans;
    FOR (i,1,300000) {
        if(adj[i].size() > 0 && freq[i] > 0)
            ans.push_back(II(freq[i], i));
    }
    sort(ans.begin(), ans.end(), greater<II> ());


    foutA << "Proportion_All(%)\tProportion_Classified(%)\tCount\tID\tName" << endl;
    for (auto i : ans) {
        if(100.0*i.fi /cntC > 0.001) {
            int u = i.se;
            int familyID = parent[u].se;
            int orderID = parent[familyID].se;
            int classID = parent[orderID].se;
            int phylumID = parent[classID].se;
            string GenusName = MName[phylumID] + " \\ " + MName[classID] + " \\ " + MName[orderID] + " \\ " + MName[familyID] + " \\ " + MName[u];
            foutA << 100.0*i.fi/cntReads << "\t" << 100.0*i.fi/cntC << "\t";
            foutA << i.fi << "\t" << u << " (G) " << GenusName << endl;


            vector<II> tmp;
            for(auto v : adj[u]) if(freq[v] > 0) {
                tmp.push_back(II(freq[v], v));
            }
            sort(tmp.begin(), tmp.end(), greater<II> ());
            for(auto j : tmp){
                foutA << 100.0*j.fi/cntReads << "\t" << 100.0*j.fi/cntC << "\t";
                foutA << j.fi << "\t" << j.se << "\t" << MName[j.se] << endl;
            }
            foutA << endl;
        }
    }
    foutA.close();

	return 0;
}

