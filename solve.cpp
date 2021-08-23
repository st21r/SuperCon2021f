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

ofstream dout("./dump.txt");

inline void dump() { dout << "\n"; }
template <typename T, typename ...U>
inline void dump(const T &t, const U &...u) {
    dout << t;
    if (sizeof...(u)) dout << " ";
    dump(u...);
}

template<typename T> inline bool chmax(T &a, T b) { if (a < b) { a = b; return 1; } return 0; }
template<typename T> inline bool chmin(T &a, T b) { if (a > b) { a = b; return 1; } return 0; }
#pragma endregion

#include <omp.h>
//#include "mpi.h"
#include "sc21.h"
#include "cssl.h"

/*
ヘッダーファイル定義一覧
#define N_GROUP 100
const double BETA=0.0002;
const double BETA2=0.000001;
const double GAMMA=0.1;
const int T=200;
int N_LINK; //辺数
int N[N_GROUP]; //各集団の人数
double I_PROB[N_GROUP]; //t日目の感染者数
int C[N_GROUP][N_GROUP]; //出力する隣接行列
double TIME0;
*/

vi randInt(int num, int r) {
    int seed = 0;
    int icon;
    double dwork[8];
    double s[num];
    vi res(num);
    c_dm_vranu5(&seed, s, num, (long)0, dwork, &icon);
    #pragma omp parallel for
    rep(i, num) {
        res[i] = (int)(s[i] * r);
        if (res[i] >= r) {
            res[i] = 0;
        }
    }
    return res;
}

vi randInt(int num, int l, int r) {
    vi res(num);
    vi &&s = randInt(num, r-l);
    #pragma omp parallel for
    rep(i, num) {
        res[i] = l + s[i];
    }
    return res;
}

void solve() {
    
}

double calcRSS() {
    double S[N_GROUP], I[N_GROUP], R[N_GROUP];
    rep(i, N_GROUP) {
        S[i] = N[i] - 1;
        I[i] = 1;
        R[i] = 0;
    }
    rep(i, N_GROUP) {
        double ni1 = BETA * S[i] * I[i], ni2 = 0, nr = GAMMA * I[i];
        rep(j, N_GROUP) {
            if (C[i][j] == 1) {
                ni2 += I[j];
            }
        }
        ni2 *= BETA2 * S[i];
        S[i] = S[i] - ni1 - ni2;
        I[i] = I[i] + ni1 + ni2 - nr;
        R[i] = R[i] + nr;
    }
    double E = 0;
    rep(i, N_GROUP) {
        double diff = I[i] - I_PROB[i];
        E += diff * diff;
    }
    return E;
}

void init() {
    rep(i, N_GROUP) {
        rep(j, N_GROUP) {
            C[i][j] = 0;
        }
    }
    int cnt = 0;
    while (cnt < N_LINK) {
        vi p = randInt(2, N_GROUP);
        int i = p[0], j = p[1];
        if (C[i][j] == 0 && i != j) {
            C[i][j] = 1;
            C[j][i] = 1;
            cnt++;
        }
    }
    dump("score:", calcRSS());
    rep(i, N_GROUP) {
        rep(j, N_GROUP) {
            dout << C[i][j];
        }
        dout << endl;
    }
}

namespace Ilag {
    

}

int main() {
    SC_input();
    init();
    SC_output();
    return 0;
}

