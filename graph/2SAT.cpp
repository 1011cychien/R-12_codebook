struct TwoSat {
  int n, N;
  vector<vector<int>> adj;
  vector<int> ans;
  TwoSat(int n) : n(n), N(n), adj(2 * n) {}
  // u == x
  void addClause(int u, bool x) { adj[2 * u + !x].push_back(2 * u + x); }
  // u == x || v == y
  void addClause(int u, bool x, int v, bool y) {
    adj[2 * u + !x].push_back(2 * v + y);
    adj[2 * v + !y].push_back(2 * u + x);
  }
  // u == x -> v == y
  void addImply(int u, bool x, int v, bool y) { addClause(u, !x, v, y); }
  void addVar() {
    adj.emplace_back(), adj.emplace_back();
    N++;
  }
  // at most one in var is true
  // adds prefix or as supplementary variables
  void atMostOne(const vector<pair<int, bool>> &vars) {
    int sz = vars.size();
    for (int i = 0; i < sz; i++) {
      addVar();
      auto [u, x] = vars[i];
      addImply(u, x, N - 1, true);
      if (i > 0) {
        addImply(N - 2, true, N - 1, true);
        addClause(u, !x, N - 2, false);
      }
    }
  }
  // does not return supplementary variables from atMostOne()
  bool satisfiable() {
    // run tarjan scc on 2 * N
    for (int i = 0; i < 2 * N; i++) { if (dfn[i] == -1) { dfs(dfs, i); }}
    for (int i = 0; i < N; i++) { if (id[2 * i] == id[2 * i + 1]) { return false; }}
    ans.resize(n);
    for (int i = 0; i < n; i++) { ans[i] = id[2 * i] > id[2 * i + 1]; }
    return true;
  }
};
