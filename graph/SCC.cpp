struct SCC {
    int n, comps = 0;
    vector<int> order, id;
    vector<vector<int>> components;
    SCC(const vector<vector<int>> &g) : n(g.size()), id(n, -1) {   
        vector<bool> used(n);
        auto dfs1 = [&](auto dfs1, int u) -> void {
            used[u] = true;
            for (int v : g[u]) {
                if (!used[v]) {
                    dfs1(dfs1, v);
                }
            }
            order.push_back(u);
        };
        for (int i = 0; i < n; ++i) {
            if (!used[i]) {
                dfs1(dfs1, i);
            }
        }
        reverse(order.begin(), order.end());
        vector<vector<int>> gr(n);
        for (int i = 0; i < n; i++) {
            for (int j : g[i]) {
                gr[j].push_back(i);
            }
        }
        used.assign(n, false);
        auto dfs2 = [&](auto dfs2, int u) -> void {
            used[u] = true;
            components.back().push_back(u);
            for (int v : gr[u]) {
                if (!used[v]) {
                    dfs2(dfs2, v);
                }
            }
        };
        for (int u : order) {
            if (!used[u]) {
                components.emplace_back();
                dfs2(dfs2, u);
                for (int v : components.back()) {
                    id[v] = comps;
                }
                comps++;
            }
        }
    }
    // the components are in topological sort order
};