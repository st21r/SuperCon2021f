vl dijkstra(const vector<vpi> &g, int s) {
    vl d(len(g), LINF);
    priority_queue<pl, vpl, greater<pl>> q;
    d[s] = 0;
    q.emplace(make_pair(0, s));
    while (!q.empty()) {
        auto [dist, v] = q.top(); q.pop();
        if (d[v] < dist) continue;
        for (auto &[nv, cost]: g[v]) {
            if (d[nv] > d[v] + cost) {
                d[nv] = d[v] + cost;
                q.emplace(make_pair(d[nv], nv));
            }
        }
    }
    return d;
}
