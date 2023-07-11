// a component contains no articulation point, so P2 is a component
// resulting bct is rooted
struct BlockCutTree {
    int n, square = 0, cur = 0;
    vector<int> low, dfn, stk;
    vector<vector<int>> g, bct;
    BlockCutTree(int n) : n(n), low(n), dfn(n, -1), g(n), bct(n) {}
    void build() { dfs(0); }
    void addEdge(int u, int v) { g[u].push_back(v), g[v].push_back(u); }
    void dfs(int u) {
        low[u] = dfn[u] = cur++;
        stk.push_back(u);
        for (auto v : g[u]) {
            if (dfn[v] == -1) {
                dfs(v);
                low[u] = min(low[u], low[v]);
                if (low[v] == dfn[u]) {
                    bct.emplace_back();
                    int x;
                    do {
                        x = stk.back();
                        stk.pop_back();
                        bct.back().push_back(x);
                    } while (x != v);
                    bct[u].push_back(n + square);
                    square++;
                }
            } else {
                low[u] = min(low[u], dfn[v]);
            }
        }
    }
};