template <typename T>
vector<T> compress(vector<T> &a) {
    vector<T> p_val = a;
    sort(all(p_val)); uni(p_val);
    for (int i = 0; i < (int)a.size(); i++) {
        a[i] = lower_bound(all(p_val), a[i]) - p_val.begin();
    }
    return p_val;
}
