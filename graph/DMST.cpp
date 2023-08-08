// DSU with rollback
template <typename Cost>
struct DMST {
    int n;
    vector<int> s, t, lc, rc, h;
    vector<Cost> c, tag;
    DMST(int n) : n(n), h(n, -1) {}
    void addEdge(int u, int v, Cost w) {
        int id = s.size();
        s.push_back(u), t.push_back(v), c.push_back(w);
        lc.push_back(-1), rc.push_back(-1);
        tag.emplace_back();
        h[v] = merge(h[v], id);
    }
    pair<Cost, vector<int>> build(int root = 0) {
        DSU d(n);
        Cost res{};
        vector<int> vis(n, -1), path(n), q(n), in(n, -1);
        vis[root] = root;
        vector<pair<int, vector<int>>> cycles;
        for (auto r = 0; r < n; ++r) {
            auto u = r, b = 0, w = -1;
            while (!~vis[u]) {
                if (!~h[u]) { return {-1, {}}; }
                push(h[u]);
                int e = h[u];
                res += c[e], tag[h[u]] -= c[e];
                h[u] = pop(h[u]);
                q[b] = e, path[b++] = u, vis[u] = r;
                u = d.find(s[e]);
                if (vis[u] == r) {
                    int cycle = -1, e = b;
                    do {
                        w = path[--b];
                        cycle = merge(cycle, h[w]);
                    } while (d.join(u, w));
                    u = d.find(u);
                    h[u] = cycle, vis[u] = -1;
                    cycles.emplace_back(u, vector<int>(q.begin() + b, q.begin() + e));
                }
            }
            for (auto i = 0; i < b; ++i) { in[d.find(t[q[i]])] = q[i]; }
        }
        reverse(cycles.begin(), cycles.end());
        for (const auto &[u, comp] : cycles) {
            int count = int(comp.size()) - 1;
            d.back(count);
            int ine = in[u];
            for (auto e : comp) { in[d.find(t[e])] = e; }
            in[d.find(t[ine])] = ine;
        }
        vector<int> par;
        par.reserve(n);
        for (auto i : in) { par.push_back(i != -1 ? s[i] : -1); }
        return {res, par};
    }
    void push(int u) {
        c[u] += tag[u];
        if (int l = lc[u]; l != -1) { tag[l] += tag[u]; }
        if (int r = rc[u]; r != -1) { tag[r] += tag[u]; }
        tag[u] = 0;
    }
    int merge(int u, int v) {
        if (u == -1 || v == -1) { return u != -1 ? u : v; }
        push(u);
        push(v);
        if (c[u] > c[v]) { swap(u, v); }
        rc[u] = merge(v, rc[u]);
        swap(lc[u], rc[u]);
        return u;
    }
    int pop(int u) {
        push(u);
        return merge(lc[u], rc[u]);
    }
};