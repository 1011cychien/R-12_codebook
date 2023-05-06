struct EBCC {
    int n;
    vector<vector<int>> g;
    vector<int> stk;
    vector<int> dfn, low, id;
    int cur, cnt;
    EBCC(const vector<vector<int>> &g, int root = 0) : n(g.size()), g(g), dfn(n, -1), low(n), id(n, -1), cur(0), cnt(0) {
        dfs(root, -1);
    }
    void dfs(int u, int p) {
        dfn[u] = low[u] = cur++;
        stk.push_back(u);
        for (auto v : g[u]) {
            if (v == p) {
                continue;
            }
            if (dfn[v] == -1) {
                dfs(v, u);
                low[u] = min(low[u], low[v]);
            } else if (id[v] == -1) {
                low[u] = min(low[u], dfn[v]);
            }
        }
        if (dfn[u] == low[u]) {
            int x;
            do {
                v = stk.back();
                id[x] = cnt;
                stk.pop_back();
            } while (x != u);
            cnt++;
        }
    }
};