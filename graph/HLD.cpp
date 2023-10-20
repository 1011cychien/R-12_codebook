struct HLD {
  int n, cur = 0;
  vector<int> sz, top, dep, par, tin, tout, seq;
  vector<vector<int>> adj;
  HLD(int n) : n(n), sz(n, 1), top(n), dep(n), par(n), tin(n), tout(n), seq(n), adj(n) {}
  void addEdge(int u, int v) { adj[u].push_back(v), adj[v].push_back(u); }
  void build(int root = 0) {
    top[root] = root, dep[root] = 0, par[root] = -1;
    dfs1(root), dfs2(root);
  }
  void dfs1(int u) {
    if (auto it = find(adj[u].begin(), adj[u].end(), par[u]); it != adj[u].end()) {
      adj[u].erase(it);
    }
    for (auto &v : adj[u]) {
      par[v] = u;
      dep[v] = dep[u] + 1;
      dfs1(v);
      sz[u] += sz[v];
      if (sz[v] > sz[adj[u][0]]) { swap(v, adj[u][0]); }
    }
  }
  void dfs2(int u) {
    tin[u] = cur++;
    seq[tin[u]] = u;
    for (auto v : adj[u]) {
      top[v] = v == adj[u][0] ? top[u] : v;
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
  int dist(int u, int v) { return dep[u] + dep[v] - 2 * dep[lca(u, v)]; }
  int jump(int u, int k) {
    if (dep[u] < k) { return -1; }
    int d = dep[u] - k;
    while (dep[top[u]] > d) { u = par[top[u]]; }
    return seq[tin[u] - dep[u] + d];
  }
  // u is v's ancestor
  bool isAncestor(int u, int v) { return tin[u] <= tin[v] && tin[v] <= tout[u]; }
  // root's parent is itself
  int rootedParent(int r, int u) {
    if (r == u) { return u; }
    if (isAncestor(r, u)) { return par[u]; }
    auto it = upper_bound(adj[u].begin(), adj[u].end(), r, [&](int x, int y) {
      return tin[x] < tin[y];
    }) - 1;
    return *it;
  }
  // rooted at u, v's subtree size
  int rootedSize(int r, int u) {
    if (r == u) { return n; }
    if (isAncestor(u, r)) { return sz[u]; }
    return n - sz[rootedParent(r, u)];
  }
  int rootedLca(int r, int a, int b) { return lca(a, b) ^ lca(a, r) ^ lca(b, r); }
};
