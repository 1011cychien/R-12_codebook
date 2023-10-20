struct BipartiteMatching {
  int n, m;
  vector<vector<int>> adj;
  vector<int> l, r, dis, cur;
  BipartiteMatching(int n, int m) : n(n), m(m), adj(n), l(n, -1), r(m, -1), dis(n), cur(n) {}
  void addEdge(int u, int v) { adj[u].push_back(v); }
  void bfs() {
    vector<int> q;
    for (int u = 0; u < n; u++) {
      if (l[u] == -1) {
        q.push_back(u), dis[u] = 0;
      } else {
        dis[u] = -1;
      }
    }
    for (int i = 0; i < int(q.size()); i++) {
      int u = q[i];
      for (auto v : adj[u]) {
        if (r[v] != -1 && dis[r[v]] == -1) {
          dis[r[v]] = dis[u] + 1;
          q.push_back(r[v]);
        }
      }
    }
  }
  bool dfs(int u) {
    for (int &i = cur[u]; i < int(adj[u].size()); i++) {
      int v = adj[u][i];
      if (r[v] == -1 || dis[r[v]] == dis[u] + 1 && dfs(r[v])) {
        l[u] = v, r[v] = u;
        return true;
      }
    }
    return false;
  }
  int maxMatching() {
    int match = 0;
    while (true) {
      bfs();
      fill(cur.begin(), cur.end(), 0);
      int cnt = 0;
      for (int u = 0; u < n; u++) {
        if (l[u] == -1) {
          cnt += dfs(u);
        }
      }
      if (cnt == 0) {
        break;
      }
      match += cnt;
    }
    return match;
  }
  auto minVertexCover() {
    vector<int> L, R;
    for (int u = 0; u < n; u++) {
      if (dis[u] == -1) {
        L.push_back(u);
      } else if (l[u] != -1) {
        R.push_back(l[u]);
      }
    }
    return pair(L, R);
  }
};
