// res : parent of each vertex in dominator tree, -1 is root, -2 if not in tree
struct DominatorTree {
  int n, cur = 0;
  vector<int> dfn, rev, fa, sdom, dom, val, rp, res;
  vector<vector<int>> adj, rdom, r;
  DominatorTree(int n) : n(n), dfn(n, -1), res(n, -2), adj(n), rdom(n), r(n) {
    rev = fa = sdom = dom = val = rp = dfn;
  }
  void addEdge(int u, int v) {
    adj[u].push_back(v);
  }
  void dfs(int u) {
    dfn[u] = cur;
    rev[cur] = u;
    fa[cur] = sdom[cur] = val[cur] = cur;
    cur++;
    for (int v : adj[u]) {
      if (dfn[v] == -1) {
        dfs(v);
        rp[dfn[v]] = dfn[u];
      }
      r[dfn[v]].push_back(dfn[u]);
    }
  }
  int find(int u, int c) {
    if (fa[u] == u) { return c != 0 ? -1 : u; }
    int p = find(fa[u], 1);
    if (p == -1) { return c != 0 ? fa[u] : val[u]; }
    if (sdom[val[u]] > sdom[val[fa[u]]]) { val[u] = val[fa[u]]; }
    fa[u] = p;
    return c != 0 ? p : val[u];
  }
  void build(int s = 0) {
    dfs(s);
    for (int i = cur - 1; i >= 0; i--) {
      for (int u : r[i]) { sdom[i] = min(sdom[i], sdom[find(u, 0)]); }
      if (i > 0) { rdom[sdom[i]].push_back(i); }
      for (int u : rdom[i]) {
        int p = find(u, 0);
        if (sdom[p] == i) {
          dom[u] = i;
        } else {
          dom[u] = p;
        }
      }
      if (i > 0) { fa[i] = rp[i]; }
    }
    res[s] = -1;
    for (int i = 1; i < cur; i++) { if (sdom[i] != dom[i]) { dom[i] = dom[dom[i]]; }}
    for (int i = 1; i < cur; i++) { res[rev[i]] = rev[dom[i]]; }
  }
};
