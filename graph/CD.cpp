vector<int> sz(n), vis(n);
auto dfs1 = [&](auto dfs1, int u, int p) -> void {
    sz[u] = 1;
    for (auto v : g[u]) {
        if (v != p && !vis[v]) {
            dfs1(dfs1, v, u);
            sz[u] += sz[v];
        }
    }
};
auto dfs2 = [&](auto dfs2, int u, int p, int tot) -> int {
    for (auto v : g[u]) {
        if (v != p && !vis[v] && 2 * sz[v] > tot) {
            return dfs2(dfs2, v, u, tot);
        }
    }
    return u;
};

auto dfs = [&](auto dfs, int cen) -> void {
    dfs1(dfs1, cen, -1);
    cen = dfs2(dfs2, cen, -1, sz[cen]);
    vis[cen] = 1;
    dfs1(dfs1, cen, -1);

    for (auto v : g[cen]) {
        if (!vis[v]) {
            dfs(dfs, v);
        }
    }
};
dfs(dfs, 0);