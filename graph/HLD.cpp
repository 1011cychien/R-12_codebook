struct HLD {
    int n, cur = 0;
    vector<int> sz, top, dep, par, tin, tout, seq;
    vector<vector<int>> g;
    HLD(int n) : n(n), sz(n, 1), top(n), dep(n), par(n), tin(n), tout(n), seq(n), g(n) {}
    void addEdge(int u, int v) { g[u].push_back(v), g[v].push_back(u); }
    void build(int root = 0) {
        top[root] = root;
        dep[root] = 0;
        par[root] = -1;
        dfs1(root);
        dfs2(root);
    }
    void dfs1(int u) {
        if (auto it = find(g[u].begin(), g[u].end(), par[u]); it != g[u].end()) {
            g[u].erase(it);
        }
        for (auto &v : g[u]) {
            par[v] = u;
            dep[v] = dep[u] + 1;
            dfs1(v);
            sz[u] += sz[v];
            if (sz[v] > sz[g[u][0]]) {
                swap(v, g[u][0]);
            }
        }
    }
    void dfs2(int u) {
        tin[u] = cur++;
        seq[tin[u]] = u;
        for (auto v : g[u]) {
            top[v] = v == g[u][0] ? top[u] : v;
            dfs2(v);
        }
        tout[u] = cur - 1;
    }
    int lca(int u, int v) {
        while (top[u] != top[v]) {
            if (dep[top[u]] > dep[top[v]]) {
                u = par[top[u]];
            } else {
                v = par[top[v]];
            }
        }
        return dep[u] < dep[v] ? u : v;
    }
    int dist(int u, int v) {
        return dep[u] + dep[v] - 2 * dep[lca(u, v)];
    }
    int jump(int u, int k) {
        if (dep[u] < k) {
            return -1;
        }
        int d = dep[u] - k;
        while (dep[top[u]] > d) {
            u = par[top[u]];
        }
        return seq[tin[u] - dep[u] + d];
    }
    // u is v's ancestor
    bool isAncestor(int u, int v) {
        return tin[u] <= tin[v] && tin[v] <= tout[u];
    }
    // rooted at u, v's parent, root's parent is itself
    int rootedParent(int u, int v) {
        if (u == v) {
            return u;
        }
        if (isAncestor(u, v)) {
            return par[v];
        }
        auto it = upper_bound(g[v].begin(), g[v].end(), u, [&](int x, int y) {
            return tin[x] < tin[y];
        }) - 1;
        return *it;
    }
    // rooted at u, v's subtree size
    int rootedSize(int u, int v) {
        if (u == v) {
            return n;
        }
        if (!isAncestor(v, u)) {
            return sz[v];
        }
        return n - sz[rootedParent(u, v)];
    }
    int rootedLca(int a, int b, int c) {
        return lca(a, b) ^ lca(b, c) ^ lca(c, a);
    }
};