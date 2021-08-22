#pragma region my_template
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <queue>
#include <string>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace std;
using ll = long long;
using pi = pair<int, int>;
using ti = tuple<int, int, int>;
using vi = vector<int>;
using vvi = vector<vi>;

#define range(i, l, r) for(int i = (int)(l); i < (int)(r); i++)
#define rrange(i, l, r) for(int i = (int)(r)-1; i >= (int)(l); i--)
#define rep(i, n) range(i, 0, n)
#define rrep(i, n) rrange(i, 0, n)
#define len(a) ((int)(a).size())
#define all(a) (a).begin(), (a).end()
#define rall(a) (a).rbegin(), (a).rend()
#define uni(a) (a).erase(unique(all(a)), (a).end());

constexpr int INF = 1e9;
constexpr ll LINF = 1e18;

template<typename T>
istream &operator >> (istream &in, vector<T> &a){
    for(T &x: a) in >> x;
    return in;
}
template<typename T, typename U>
istream &operator >> (istream &in, pair<T, U> &a){
    in >> a.first >> a.second;
    return in;
}
template<typename T>
ostream &operator << (ostream &out, const vector<T> &a) {
    rep(i, len(a)) out << a[i] << (i == len(a)-1 ? "" : " ");
    return out;
}
template<typename T, typename U>
ostream &operator << (ostream &out, const pair<T, U> &a){
    out << a.first << " " << a.second;
    return out;
}
inline void print() { cout << "\n"; }
template <typename T, typename ...U>
inline void print(const T &t, const U &...u) {
    cout << t;
    if (sizeof...(u)) cout << " ";
    print(u...);
}
#ifdef LOCAL
ofstream dout("./dump.txt");
#else
ofstream dout("/dev/null");
#endif
inline void dump() { dout << "\n"; }
template <typename T, typename ...U>
inline void dump(const T &t, const U &...u) {
    dout << t;
    if (sizeof...(u)) dout << " ";
    dump(u...);
}

template<typename T> inline bool chmax(T &a, T b) { if (a < b) { a = b; return 1; } return 0; }
template<typename T> inline bool chmin(T &a, T b) { if (a > b) { a = b; return 1; } return 0; }

namespace util {
    struct Timer {
        static const uint64_t CYCLES_PER_SEC = 3e9;
        uint64_t start;
    
        Timer() : start{} { reset(); }
    
        void reset() { start = get_cycle(); }
    
        inline double get() const { return (double) get_cycle() / CYCLES_PER_SEC; }

    private:
        inline uint64_t get_cycle() const {
            unsigned low, high;
            __asm__ volatile ("rdtsc" : "=a" (low), "=d" (high));
            return (((uint64_t) low) | ((uint64_t) high << 32ull)) - start;
        }
    };

    class XorShift {
        unsigned x, y, z, w; 
    public:    
        XorShift() {
            x = 123456789;
            y = 362436069;
            z = 521288629;
            w = 88675123;
        }
        inline unsigned next() {
            unsigned t = x ^ (x << 11);
            x = y; y = z; z = w;
            return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
        }
        unsigned nextInt(int n) {
            return next() % n;
        }
        unsigned nextInt(int l, int r) {
            return l + (next() % (r - l));
        }
        double nextDouble() {
            return next() / 4294967295.0;
        }
    };
}

#pragma endregion

void init() {
    
}

int main() {
    cin.tie(0);
    ios::sync_with_stdio(false);
    init();
    return 0;
}
