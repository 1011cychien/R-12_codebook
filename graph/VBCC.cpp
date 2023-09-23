// is articulation point if appear in >= 2 comps
auto dfs = [&](auto dfs, int u, int p) -> void {
    dfn[u] = low[u] = T++;
    for (auto v : adj[u]) {
        if (v == p) { continue; }
        if (dfn[v] == -1) {
            stk.push_back(v);
            dfs(dfs, v, u);
            low[u] = min(low[u], low[v]);
            if (low[v] >= dfn[u]) {
                comps.emplace_back();
                int x;
                do {
                    x = stk.back();
                    cnt[x]++;
                    stk.pop_back();
                } while (x != v);
                comps.back().push_back(u);
                cnt[u]++;
            }
        } else {
            low[u] = min(low[u], dfn[v]);
        }
    }
};
for (int i = 0; i < n; i++) {
    if (!adj[i].empty()) {
        dfs(dfs, i, -1);
    } else {
        comps.push_back({i});
    }
}