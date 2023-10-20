// O(V ^ 3)
template <typename F>
struct GlobalMinCut {
  static constexpr int INF = numeric_limits<F>::max() / 2;
  int n;
  vector<int> vis, wei;
  vector<vector<int>> adj;
  GlobalMinCut(int n) : n(n), vis(n), wei(n), adj(n, vector<int>(n)) {}
  void addEdge(int u, int v, int w){
    adj[u][v] += w;
    adj[v][u] += w;
  }
  int solve() {
    int sz = n;
    int res = INF, x = -1, y = -1;
    auto search = [&]() {
      fill(vis.begin(), vis.begin() + sz, 0);
      fill(wei.begin(), wei.begin() + sz, 0);
      x = y = -1;
      int mx, cur;
      for (int i = 0; i < sz; i++) {
        mx = -1, cur = 0;
        for (int j = 0; j < sz; j++) {
          if (wei[j] > mx) {
            mx = wei[j], cur = j;
          }
        }
        vis[cur] = 1, wei[cur] = -1;
        x = y;
        y = cur;
        for (int j = 0; j < sz; j++) {
          if (!vis[j]) {
            wei[j] += adj[cur][j];
          }
        }
      }
      return mx;
    };
    while (sz > 1) {
      res = min(res, search());
      for (int i = 0; i < sz; i++) {
        adj[x][i] += adj[y][i];
        adj[i][x] = adj[x][i];
      }
      for (int i = 0; i < sz; i++) {
        adj[y][i] = adj[sz - 1][i];
        adj[i][y] = adj[i][sz - 1];
      }
      sz--;
    }
    return res;
  }
};
