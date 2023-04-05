vector<int> tin(n), tout(n), sz(n, 1), dep(n), top(n), par(n);
// root = 0 ? 
top[0] = 0, par[0] = -1;
auto dfs1 = [&](auto dfs1, int u, int p) -> void {
    if (p != -1) {
        g[u].erase(find(g[u].begin(), g[u].end(), p));
    }
    for (auto &v : g[u]) {
        par[v] = u;
        dep[v] = dep[u] + 1;
        dfs1(dfs1, v, u);
        sz[u] += sz[v];
        if (sz[v] > sz[g[u][0]]) {
            swap(v, g[u][0]);
        }
    }
};
dfs1(dfs1, 0, -1);
int T = 0;
auto dfs2 = [&](auto dfs2, int u) -> void {
    tin[u] = T++;
    for (auto v : g[u]) {
        top[v] = v == g[u][0] ? top[u] : v;
        dfs2(dfs2, v);
    }
    tout[u] = T - 1;
};
dfs2(dfs2, 0);
auto lca = [&](int u, int v) {
    while (top[u] != top[v]) {
        if (dep[top[u]] < dep[top[v]]) {
            swap(u, v);
        }
        u = par[top[u]];
    }
    return dep[u] < dep[v] ? u : v;
};