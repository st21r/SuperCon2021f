template<typename T>
class BIT {
    int N;
    vector<T> data;
public:
    BIT(int n) : N(n+1), data(N) {}
    T sum(int r) {
        T res = 0;
        for (int i = r; i > 0; i -= i & -i) {
            res += data[i];
        }
        return res;
    }
    T sum() { return sum(N-1); }
    T sum(int l, int r) {
        return sum(r) - sum(l);
    }
    void add(int idx, T x = 1) {
        for(int i = idx+1; i < N; i += i & -i) {
            data[i] += x;
        }
    }
    int lower_bound(T val) {
        int x = 0, r = 1;
        while (r < N) r = r << 1;
        for (int k = r; k > 0; k = k >> 1) {
            if (x + k < N && data[x + k] < val) {
                val -= data[x + k];
                x += k;
            }
        }
        return x;
    }
};
