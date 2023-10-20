// DSU
struct ETCC {
  int n, cnt = 0;
  vector<vector<int>> adj, comps;
  vector<int> in, out, low, up, nx, id;
  ETCC(int n) : n(n), adj(n), in(n, -1), out(in), low(n), up(n), nx(in), id(in) {}
  void addEdge(int u, int v) {
    adj[u].push_back(v);
    adj[v].push_back(u);
  }
  void build() {
    int T = 0;
    DSU d(n);
    auto merge = [&](int u, int v) {
      d.join(u, v);
      up[u] += up[v];
    };
    auto dfs = [&](auto dfs, int u, int p) -> void {
      in[u] = low[u] = T++;
      for (auto v : adj[u]) {
        if (v == u) { continue; }
        if (v == p) {
          p = -1;
          continue;
        }
        if (in[v] == -1) {
          dfs(dfs, v, u);
          if (nx[v] == -1 && up[v] <= 1) {
            up[u] += up[v];
            low[u] = min(low[u], low[v]);
            continue;
          }
          if (up[v] == 0) { v = nx[v]; }
          if (low[u] > low[v]) { low[u] = low[v], swap(nx[u], v); }
          while (v != -1) { merge(u, v); v = nx[v]; }
        } else if (in[v] < in[u]) {
          low[u] = min(low[u], in[v]);
          up[u]++;
        } else {
          for (int &x = nx[u]; x != -1 && in[x] <= in[v] && in[v] < out[x]; x = nx[x]) {
            merge(u, x);
          }
          up[u]--;
        }
      }
      out[u] = T;
    };
    for (int i = 0; i < n; i++) { if (in[i] == -1) { dfs(dfs, i, -1); }}
    for (int i = 0; i < n; i++) { if (d.find(i) == i) { id[i] = cnt++; }}
    comps.resize(cnt);
    for (int i = 0; i < n; i++) { comps[id[d.find(i)]].push_back(i); }
  }
};
