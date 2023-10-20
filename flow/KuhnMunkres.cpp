// need perfect matching or not : w intialize with -INF / 0
template <typename Cost> 
struct KM {
  static constexpr Cost INF = numeric_limits<Cost>::max() / 2;
  int n;
  vector<Cost> hl, hr, slk;
  vector<int> l, r, pre, vl, vr;
  queue<int> q;
  vector<vector<Cost>> w;
  KM(int n) : n(n), hl(n), hr(n), slk(n), l(n, -1), r(n, -1), pre(n), vl(n), vr(n),
    w(n, vector<Cost>(n, -INF)) {}
  bool check(int x) {
    vl[x] = true;
    if (l[x] != -1) {
      q.push(l[x]);
      return vr[l[x]] = true;
    }
    while (x != -1) { swap(x, r[l[x] = pre[x]]); }
    return false;
  }
  void bfs(int s) {
    fill(slk.begin(), slk.end(), INF);
    fill(vl.begin(), vl.end(), false);
    fill(vr.begin(), vr.end(), false);
    q = {};
    q.push(s);
    vr[s] = true;
    while (true) {
      Cost d;
      while (!q.empty()) {
        int y = q.front();
        q.pop();
        for (int x = 0; x < n; ++x) {
          if (!vl[x] && slk[x] >= (d = hl[x] + hr[y] - w[x][y])) {
            pre[x] = y;
            if (d != 0) {
              slk[x] = d;
            } else if (!check(x)) {
              return;
            }
          }
        }
      }
      d = INF;
      for (int x = 0; x < n; ++x) { if (!vl[x] && d > slk[x]) { d = slk[x]; }}
      for (int x = 0; x < n; ++x) {
        if (vl[x]) {
          hl[x] += d;
        } else {
          slk[x] -= d;
        }
        if (vr[x]) { hr[x] -= d; }
      }
      for (int x = 0; x < n; ++x) { if (!vl[x] && !slk[x] && !check(x)) { return; }}
    }
  }
  void addEdge(int u, int v, Cost x) { w[u][v] = max(w[u][v], x); }
  Cost solve() {
    for (int i = 0; i < n; ++i) { hl[i] = *max_element(w[i].begin(), w[i].end()); }
    for (int i = 0; i < n; ++i) { bfs(i); }
    Cost res = 0;
    for (int i = 0; i < n; ++i) { res += w[i][l[i]]; }
    return res;
  }
};
