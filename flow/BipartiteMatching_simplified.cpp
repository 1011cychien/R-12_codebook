struct BipartiteMatching {
    int n, m;
    vector<vector<int>> g;
    vector<int> l, r, dis, cur;
    BipartiteMatching(int n, int m) : n(n), m(m), g(n), l(n, -1), r(m, -1), dis(n), cur(n) {}
    // come on, you know how to write this
    void addEdge(int u, int v) {}
    void bfs() {}
    bool dfs(int u) {}
    int maxMatching() {}
    auto minVertexCover() {
        vector<int> L, R;
        for (int u = 0; u < n; u++) {
            if (dis[u] == -1) {
                L.push_back(u);
            } else if (l[u] != -1) {
                R.push_back(l[u]);
            }
        }
        return pair(L, R);
    }
};