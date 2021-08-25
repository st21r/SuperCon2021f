#pragma region my_template
#include <algorithm>
#include <cmath>
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
inline void print() { cout << endl; }
template <typename T, typename ...U>
inline void print(const T &t, const U &...u) {
    cout << t;
    if (sizeof...(u)) cout << " ";
    print(u...);
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

class XorShift {
    unsigned x, y, z, w; 
public:    
    XorShift(unsigned seed) {
        x = 123456789 + seed;
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
#pragma endregion

#include <omp.h>
#include "mpi.h"
#include "cssl.h"
#include "sc21.h"

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

inline double get_time() {
    return omp_get_wtime() - TIME0;
}

XorShift rnd = XorShift(0);

int id, procs;

constexpr int M = N_GROUP*(N_GROUP-1)/2;
vector<pi> idx_to_edge(M);
double score;
ll step_num = 0;

void debug_output() {
    print("id:", id, "step:", step_num, "score:", score, "time:", get_time());
}

double eval(bool output=false) {
    vector<double> S(N_GROUP);
    rep(i, N_GROUP) {
        S[i] = (double)N[i];
    }
    vector<double> I(N_GROUP, 0.0);
    I[0]++; S[0]--;
    
    rep(_, T) {
        vector<double> PS = S, PI = I;
        rep(i, N_GROUP) {
            double ni1 = BETA * PS[i] * PI[i], ni2 = 0;
            double nr = GAMMA * PI[i];
            rep(j, N_GROUP) {
                ni2 += PI[j] * C[i][j];
            }
            ni2 *= BETA2 * PS[i];
            S[i] = PS[i] - ni1 - ni2;
            I[i] = PI[i] + ni1 + ni2 - nr;
        }
    }
    double E = 0;
    rep(i, N_GROUP) {
        double diff = I[i] - I_PROB[i];
        if (output) print(i, I[i], I_PROB[i], diff);
        E += diff * diff;
    }
    return E;
}



constexpr double TIME_LIM = 20.0;
constexpr double OUTPUT_CYCLE = 5.0;

constexpr double START_TEMP = 1000;
constexpr double END_TEMP = 0;

void gather() {
    int send = score;
    int *res;
    res = (int *)malloc(procs*sizeof(int));
    MPI_Gather(&send, 1, MPI_INT, res, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int min_id = -1;
    if (id == 0) {
        int min_score = 1e9;
        rep(i, procs) {
            if (res[i] < min_score) {
                min_score = res[i];
                min_id = i;
            }
        }
    }
    MPI_Bcast(&min_id, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (id == min_id) {
        debug_output();
    }
}

void solve() {
    vi output_flags((int)(TIME_LIM/OUTPUT_CYCLE) + 1, 0);
    while (true) {
        double cur_time = get_time();
        if (cur_time > TIME_LIM) break;
        double temp = START_TEMP + (END_TEMP - START_TEMP) * cur_time / TIME_LIM;
        
        int fi, fj, ti, tj;
        if (rnd.nextDouble() < 0.5) {
            //遷移1
            while (true) {
                tie(fi, fj) = idx_to_edge[rnd.nextInt(M)];
                tie(ti, tj) = idx_to_edge[rnd.nextInt(M)];
                if (C[fi][fj] != C[ti][tj]) {
                    if (C[fi][fj] == 0) {
                        swap(fi, ti); swap(fj, tj);
                    }
                    break;
                }
            }
        } else {
            //遷移2
            while (true) {
                int i = rnd.nextInt(N_GROUP);
                vi a0, a1;
                rep(j, N_GROUP) {
                    if (i == j) continue;
                    if (C[i][j] == 0) {
                        a0.push_back(j);
                    } else {
                        a1.push_back(j);
                    }
                }
                if (!a0.empty() && !a1.empty()) {
                    fi = i, ti = i;
                    fj = a1[rnd.nextInt(len(a1))];
                    tj = a0[rnd.nextInt(len(a0))];
                    break;
                }
            }
        }
        swap(C[fi][fj], C[ti][tj]);
        swap(C[fj][fi], C[tj][ti]);
        double nscore = eval();
        double prob = exp((score - nscore) / temp);
        if (rnd.nextDouble() < prob) {
            score = nscore;
        } else {
            swap(C[fi][fj], C[ti][tj]);
            swap(C[fj][fi], C[tj][ti]);
        }
        int o_idx = (int)(cur_time / OUTPUT_CYCLE);
        if (output_flags[o_idx] == 0) {
            gather();
            output_flags[o_idx] = 1;
        }
        step_num++;
    }
}

bool check_ans() {
    int cnt = 0;
    rep(i, N_GROUP) {
        range(j, i+1, N_GROUP) {
            if (C[i][j] != C[j][i]) return false;
            if (C[i][j]) cnt++;
        }
    }
    return cnt == N_LINK;
}

void init() {
    int cnt = 0;
    rep(i, N_GROUP) {
        rep(j, N_GROUP) {
            C[i][j] = 0;
            if (i < j) {
                idx_to_edge[cnt] = {i, j};
                cnt++;
            }
        }
    }
    //初期解の生成
    cnt = 0;
    while (cnt < N_LINK) {
        int r = rnd.nextInt(M);
        int i, j; tie(i, j) = idx_to_edge[r];
        if (C[i][j] == 0) {
            C[i][j] = 1;
            C[j][i] = 1;
            cnt++;
        }
    }
    score = eval();
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    SC_input();
    MPI_Comm_size(MPI_COMM_WORLD, &procs); 
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    rnd = XorShift(id);
    
    init();
    solve();
    
    gather();
    
    MPI_Finalize();
    //SC_output();
    return 0;
}

