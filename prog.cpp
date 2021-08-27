#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>

using namespace std;
using ll = long long;
using pi = pair<int, int>;
using ti = tuple<int, int, int>;
using vi = vector<int>;
using vvi = vector<vi>;

#define rrange(i, l, r) for(int i = (int)(r)-1; i >= (int)(l); i--)
#define rep(i, n) range(i, 0, n)
#define rrep(i, n) rrange(i, 0, n)
#define len(a) ((int)(a).size())
#define all(a) (a).begin(), (a).end()
#define rall(a) (a).rbegin(), (a).rend()

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

constexpr double TIME_LIM = 299.8;

inline double get_time() {
    return omp_get_wtime() - TIME0;
}

XorShift rnd = XorShift(0);

int id, procs;

constexpr int M = N_GROUP*(N_GROUP-1)/2;
vector<pi> idx_to_edge(M);
double score;
ll step_num = 0;

double temp;
double bscore_all;
double bid_all;

double eval() {
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
        E += diff * diff;
    }
    return E;
}

void gather() {
    int cur_time = get_time();
    double send = score;
    double *res;
    res = (double *)malloc(procs*sizeof(double));
    MPI_Gather(&send, 1, MPI_DOUBLE, res, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    double best_score = 1e9;
    double send_score = 1e9;
    int best_id = -1;
    int send_id = -1;
    
    if (id == 0) {
        vector<pair<double, int>> scores(procs);
        rep(i, procs) {
            scores[i] = {res[i], i};
        }
        sort(all(scores));
        best_score = scores[0].first;
        best_id = scores[0].second;
        // スコアの悪いとこに渡す候補数
        const int top_id_num = max(16 * (int)(1 - cur_time / TIME_LIM), 1);
        int sid = rnd.nextInt(top_id_num);
        send_score = scores[sid].first;
        send_id = scores[sid].second;
    }
    MPI_Bcast(&best_score, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&send_score, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&best_id, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&send_id, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    int C_flat[N_GROUP*N_GROUP];
    
    if (id == best_id && cur_time >= TIME_LIM-2.0) {
        SC_output();
    }
    if (id == send_id) {
        int idx = 0;
        rep(i, N_GROUP) rep(j, N_GROUP) {
            C_flat[idx] = C[i][j];
            idx++;
        }
    }
    MPI_Bcast(&C_flat[0], N_GROUP*N_GROUP, MPI_INT, send_id, MPI_COMM_WORLD);
    // 送られてきたスコアより自分のスコアが1.2倍以上悪かった時に確率で受け取る
    double s = score / send_score - 1.2;
    if (rnd.nextDouble() < s * 0.1 * (cur_time / TIME_LIM)) {
        int idx = 0;
        rep(i, N_GROUP) rep(j, N_GROUP) {
            C[i][j] = C_flat[idx];
            idx++;
        }
        score = send_score;
    }
    
    
    bscore_all = best_score;
    bid_all = best_id;
}

constexpr double GATHER_CYCLE = 1.0;

constexpr double START_TEMP = 75.0;

void solve() {
    vi output_flags((int)(TIME_LIM/GATHER_CYCLE) + 1, 0);
    while (true) {
        double cur_time = get_time();
        if (cur_time > TIME_LIM) break;
        double s = min(pow(bscore_all / score, 0.5), 1.0);
        if (id == bid_all) {
            s = 1;
        }
        double t = (cur_time / TIME_LIM);
        double f = pow(t, 0.5);
        temp = (START_TEMP - START_TEMP * f);
        
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
        int o_idx = (int)(cur_time / GATHER_CYCLE);
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
    bscore_all = score / 10;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    SC_input();
    MPI_Comm_size(MPI_COMM_WORLD, &procs); 
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    bid_all = id;
    rnd = XorShift(id);
    
    init();
    solve();
    
    gather();
    if (id == bid_all) {
        eval();
        SC_output();
    }
    
    MPI_Finalize();
    return 0;
}
