#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <sys/time.h>
#include <sys/resource.h>
#define DEBUG(a) {cerr << #a << ": " << (a) << endl; fflush(stderr); }

using namespace std;

template<class T> int getbit(T s, int i) { return (s >> i) & 1; }
template<class T> T onbit(T s, int i) { return s | (T(1) << i); }
template<class T> T offbit(T s, int i) { return s & (~(T(1) << i)); }
template<class T> int cntbit(T s) { return __builtin_popcount(s);}

#define LL long long
#define ULL unsigned long long
typedef pair<int, int> II;
typedef pair<uint32_t, uint32_t> II32;
typedef pair<II, int> III;
typedef pair<string, int> IIS;
typedef pair<string, II> IISS;

inline void printRam() {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    std::cerr << "Max ram (in kilobytes): " << ru.ru_maxrss << std::endl;
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

    void print_time(std::string s) {
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
    for(size_t i = 0; i < s.size(); ++i)
        ans = ans*10 + s[i] - 48;
    return ans;
}

uint64_t toNumDNA(string &s, int a, int len) {
    uint64_t ans = 0;
    for (size_t i = a; i < a+len; ++i) {
        if(s[i] == 'x') return 0;
        ans <<= 2;
        ans |= get_code(s[i]);
    }
    return ans;
}

uint64_t reverseMask(uint64_t _ikmer, int m_k) {
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
