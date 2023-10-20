// edu13F MLE with non-deleted pointers
// [) interval because of negative numbers
constexpr i64 INF64 = 4e18;
struct Line {
  i64 a = -INF64, b = -INF64;
  i64 operator()(i64 x) const {
    if (a == -INF64 && b == -INF64) {
      return -INF64;
    } else {
      return a * x + b;
    }
  }
};
constexpr int INF32 = 1e9;
struct LiChao {
  static constexpr int N = 5e6;
  array<Line, N> st;
  array<int, N> lc, rc;
  int n = 0;
  void clear() { n = 0; node(); }
  int node() {
    st[n] = {};
    lc[n] = rc[n] = -1;
    return n++;
  }
  void add(int id, int l, int r, Line line) {
    int m = (l + r) / 2;
    bool lcp = st[id](l) < line(l);
    bool mcp = st[id](m) < line(m);
    if (mcp) { swap(st[id], line); }
    if (r - l == 1) { return; }
    if (lcp != mcp) {
      if (lc[id] == -1) {
        lc[id] = node();
      }
      add(lc[id], l, m, line);
    } else {
      if (rc[id] == -1) {
        rc[id] = node();
      }
      add(rc[id], m, r, line);
    }
  }
  void add(Line line, int l = -INF32 - 1, int r = INF32 + 1) {
    add(0, l, r, line);
  }
  i64 query(int id, int l, int r, i64 x) {
    i64 res = st[id](x);
    if (r - l == 1) { return res; }
    int m = (l + r) / 2;
    if (x < m && lc[id] != -1) {
      res = max(res, query(lc[id], l, m, x));
    } else if (x >= m && rc[id] != -1) {
      res = max(res, query(rc[id], m, r, x));
    }
    return res;
  }
  i64 query(i64 x, int l = -INF32 - 1, int r = INF32 + 1) {
    return query(0, l, r, x);
  }
};
