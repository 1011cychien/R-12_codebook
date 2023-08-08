auto gomory(int n, vector<array<int, 3>> e) {
    Flow<int, int> mf(n);
    for (auto [u, v, c] : e) { mf.addEdge(u, v, c, c); }
    vector<array<int, 3>> res;
    vector<int> p(n);
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < int(e.size()); j++) { mf.e[j << 1].cap = mf.e[j << 1 | 1].cap = e[j][2]; }
        int f = mf.maxFlow(i, p[i]);
        auto cut = mf.minCut();
        for (int j = i + 1; j < n; j++) { if (cut[i] == cut[j] && p[i] == p[j]) { p[j] = i; }}
        res.push_back({f, i, p[i]});
    }
    return res;
}