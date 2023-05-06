struct HLD {
    int n;
    vector<int> sz, top, dep, par, tin, tout, seq;
    vector<vector<int>> g;
    int cur;
    HLD(const vector<vector<int>> &g, int root = 0) : n(g.size()), sz(n, 1), top(n), dep(n), par(n), tin(n), tout(n), seq(n), cur(0), g(g) {
        top[root] = root;
        dep[root] = 0;
        par[root] = -1;
        dfs1(root);
        dfs2(root);
    }
    void dfs1(int u) {
        if (par[u] != -1) {
            g[u].erase(find(g[u].begin(), g[u].end(), par[u]));
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
        tout[u] = cur;
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
        return tin[u] <= tin[v] && tin[v] < tout[u];
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