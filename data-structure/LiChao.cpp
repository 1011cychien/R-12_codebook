constexpr i64 INF = 4e18;
struct Line {
  i64 a, b;
  Line() : a(0), b(INF) {}
  Line(i64 a, i64 b) : a(a), b(b) {}
  i64 operator()(i64 x) { return a * x + b; }
};
// [, ) !!!!!!!!!!!
struct Lichao {
  int n;
  vector<int> vals;
  vector<Line> lines;
  Lichao() {}
  void init(const vector<int> &v) {
    n = v.size();
    vals = v;
    sort(vals.begin(), vals.end());
    vals.erase(unique(vals.begin(), vals.end()), vals.end());
    lines.assign(4 * n, {});
  }
  int get(int x) { return lower_bound(vals.begin(), vals.end(), x) - vals.begin(); }
  void apply(Line p, int id, int l, int r) {
    Line &q = lines[id];
    if (p(vals[l]) < q(vals[l])) { swap(p, q); }
    if (l + 1 == r) { return; }
    int m = l + r >> 1;
    if (p(vals[m]) < q(vals[m])) {
      swap(p, q);
      apply(p, id << 1, l, m);
    } else {
      apply(p, id << 1 | 1, m, r);
    }
  }
  void add(int ql, int qr, Line p) {
    ql = get(ql), qr = get(qr);
    auto go = [&](auto go, int id, int l, int r) -> void {
      if (qr <= l || r <= ql) { return; }
      if (ql <= l && r <= qr) {
        apply(p, id, l, r);
        return;
      }
      int m = l + r >> 1;
      go(go, id << 1, l, m);
      go(go, id << 1 | 1, m, r);
    };
    go(go, 1, 0, n);
  }
  i64 query(int p) {
    p = get(p);
    auto go = [&](auto go, int id, int l, int r) -> i64 {
      if (l + 1 == r) { return lines[id](vals[p]); }
      int m = l + r >> 1;
      return min(lines[id](vals[p]), p < m ? go(go, id << 1, l, m) : go(go, id << 1 | 1, m, r));
    };
    return go(go, 1, 0, n);
  }
};
