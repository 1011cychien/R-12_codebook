// block cut tree, directed
// note that edge to father is not ignored in this implementation! no problem!
int square = 0;
vector<int> low(n), dfn(n, -1), stk;
vector<vector<int>> bct(n);
auto tarjan = [&](auto tarjan, int u) -> void {
    static int T = 0;
    dfn[u] = low[u] = T++;
    stk.push_back(u);
    for (auto v : g[u]) {
        if (dfn[v] == -1) {
            tarjan(tarjan, v);
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
};
tarjan(tarjan, 0);