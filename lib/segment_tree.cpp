class SegmentTree {
    const ll IDE = LINF;
    inline ll f(ll x, ll y) {
        return min(x, y);
    }
    
    int n;
    vl data;
public:
    SegmentTree(int size) {
        n = 1;
        while (n < size) n <<= 1;
        data.assign(2*n, IDE);
    }

    void build(vl a) {
        rep(i, len(data)) {
            data[n+i] = a[i];
        }
        rrep(i, n) {
            data[i] = f(data[i*2], data[i*2+1]);
        }
    }
    void update(int i, ll x) {
        i += n;
        data[i] = x;
        while (i > 1) {
            i >>= 1;
            data[i] = f(data[i*2], data[i*2+1]);
        }
    }
    ll query(int l, int r) {
        l += n; r += n;
        ll res = IDE;
        while (l < r) {
            if (l & 1) {
                res = f(res, data[l]);
                l += 1;
            }
            if (r & 1) {
                res = f(res, data[r-1]);
            }
            l >>= 1; r >>= 1;
        }
        return res;
    }
    ll get(int i) {
        return data[i+n];
    }
};
