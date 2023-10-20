struct SCC {
  int n, cnt = 0, cur = 0;
  vector<int> id, dfn, low, stk;
  vector<vector<int>> adj, comps;
  void addEdge(int u, int v) { adj[u].push_back(v); }
  SCC(int n) : n(n), id(n, -1), dfn(n, -1), low(n, -1), adj(n) {}
  void build() {
    auto dfs = [&](auto dfs, int u) -> void {
      dfn[u] = low[u] = cur++;
      stk.push_back(u);
      for (auto v : adj[u]) {
        if (dfn[v] == -1) {
          dfs(dfs, v);
          low[u] = min(low[u], low[v]);
        } else if (id[v] == -1) {
          low[u] = min(low[u], dfn[v]);
        }
      }
      if (dfn[u] == low[u]) {
        int v;
        comps.emplace_back();
        do {
          v = stk.back();
          comps.back().push_back(v);
          id[v] = cnt;
          stk.pop_back();
        } while (u != v);
        cnt++;
      }
    };
    for (int i = 0; i < n; i++) { if (dfn[i] == -1) { dfs(dfs, i); }}
    for (int i = 0; i < n; i++) { id[i] = cnt - 1 - id[i]; }
    reverse(comps.begin(), comps.end());
  }
  // the comps are in topological sorted order
};
