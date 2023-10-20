struct Line {
  // kx + b
  mutable i64 k, b, p;
  bool operator<(const Line& o) const { return k < o.k; }
  bool operator<(i64 x) const { return p < x; }
};
struct DynamicConvexHullMax : multiset<Line, less<>> {
  // (for doubles, use INF = 1/.0, div(a,b) = a/b)
  static constexpr i64 INF = numeric_limits<i64>::max();
  i64 div(i64 a, i64 b) {
    // floor
    return a / b - ((a ^ b) < 0 && a % b);
  }
  bool isect(iterator x, iterator y) {
    if (y == end()) return x->p = INF, 0;
    if (x->k == y->k) x->p = x->b > y->b ? INF : -INF;
    else x->p = div(y->b - x->b, x->k - y->k);
    return x->p >= y->p;
  }
  void add(i64 k, i64 b) {
    auto z = insert({k, b, 0}), y = z++, x = y;
    while (isect(y, z)) z = erase(z);
    if (x != begin() && isect(--x, y)) isect(x, y = erase(y));
    while ((y = x) != begin() && (--x)->p >= y->p)
      isect(x, erase(y));
  }
  i64 query(i64 x) {
    if (empty()) {
      return -INF;
    }
    auto l = *lower_bound(x);
    return l.k * x + l.b;
  }
};
