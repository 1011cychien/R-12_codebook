struct SegmentTree {
  int n;
  struct node {
    i64 mx1, mx2, mxc;
    i64 mn1, mn2, mnc;
    i64 add;
    i64 sum;
    node(i64 v = 0) {
      mx1 = mn1 = sum = v;
      mxc = mnc = 1;
      add = 0;
      mx2 = -9e18, mn2 = 9e18;
    }
  };
  vector<node> t;
  // build
  void push(int id, int l, int r) {
    auto& c = t[id];
    int m = l + r >> 1;
    if (c.add != 0) {
      apply_add(id << 1, l, m, c.add);
      apply_add(id << 1 | 1, m + 1, r, c.add);
      c.add = 0;
    }
    apply_min(id << 1, l, m, c.mn1);
    apply_min(id << 1 | 1, m + 1, r, c.mn1);
    apply_max(id << 1, l, m, c.mx1);
    apply_max(id << 1 | 1, m + 1, r, c.mx1);
  }
  void apply_add(int id, int l, int r, i64 v) {
    if (v == 0) {
      return;
    }
    auto& c = t[id];
    c.add += v;
    c.sum += v * (r - l + 1);
    c.mx1 += v;
    c.mn1 += v;
    if (c.mx2 != -9e18) {
      c.mx2 += v;
    }
    if (c.mn2 != 9e18) {
      c.mn2 += v;
    }
  }
  void apply_min(int id, int l, int r, i64 v) {
    auto& c = t[id];
    if (v <= c.mn1) {
      return;
    }
    c.sum -= c.mn1 * c.mnc;
    c.mn1 = v;
    c.sum += c.mn1 * c.mnc;
    if (l == r || v >= c.mx1) {
      c.mx1 = v;
    } else if (v > c.mx2) {
      c.mx2 = v;
    }
  }
  void apply_max(int id, int l, int r, i64 v) {
    auto& c = t[id];
    if (v >= c.mx1) {
      return;
    }
    c.sum -= c.mx1 * c.mxc;
    c.mx1 = v;
    c.sum += c.mx1 * c.mxc;
    if (l == r || v <= c.mn1) {
      c.mn1 = v;
    } else if (v < c.mn2) {
      c.mn2 = v;
    }
  }
  void pull(int id) {
    auto &c = t[id], &lc = t[id << 1], &rc = t[id << 1 | 1];
    c.sum = lc.sum + rc.sum;
    if (lc.mn1 == rc.mn1) {
      c.mn1 = lc.mn1;
      c.mn2 = min(lc.mn2, rc.mn2);
      c.mnc = lc.mnc + rc.mnc;
    } else if (lc.mn1 < rc.mn1) {
      c.mn1 = lc.mn1;
      c.mn2 = min(lc.mn2, rc.mn1);
      c.mnc = lc.mnc;
    } else {
      c.mn1 = rc.mn1;
      c.mn2 = min(lc.mn1, rc.mn2);
      c.mnc = rc.mnc;
    }
    if (lc.mx1 == rc.mx1) {
      c.mx1 = lc.mx1;
      c.mx2 = max(lc.mx2, rc.mx2);
      c.mxc = lc.mxc + rc.mxc;
    } else if (lc.mx1 > rc.mx1) {
      c.mx1 = lc.mx1;
      c.mx2 = max(lc.mx2, rc.mx1);
      c.mxc = lc.mxc;
    } else {
      c.mx1 = rc.mx1;
      c.mx2 = max(lc.mx1, rc.mx2);
      c.mxc = rc.mxc;
    }
  }
  void range_chmin(int id, int l, int r, int ql, int qr, i64 v) {
    if (r < ql || l > qr || v >= t[id].mx1) {
      return;
    }
    if (ql <= l && r <= qr && v > t[id].mx2) {
      apply_max(id, l, r, v);
      return;
    }
    push(id, l, r);
    int m = l + r >> 1;
    range_chmin(id << 1, l, m, ql, qr, v);
    range_chmin(id << 1 | 1, m + 1, r, ql, qr, v);
    pull(id);
  }
  void range_chmin(int ql, int qr, i64 v) {
    range_chmin(1, 0, n - 1, ql, qr, v);
  }
  void range_chmax(int id, int l, int r, int ql, int qr, i64 v) {
    if (r < ql || l > qr || v <= t[id].mn1) {
      return;
    }
    if (ql <= l && r <= qr && v < t[id].mn2) {
      apply_min(id, l, r, v);
      return;
    }
    push(id, l, r);
    int m = l + r >> 1;
    range_chmax(id << 1, l, m, ql, qr, v);
    range_chmax(id << 1 | 1, m + 1, r, ql, qr, v);
    pull(id);
  }
  void range_chmax(int ql, int qr, i64 v) {
    range_chmax(1, 0, n - 1, ql, qr, v);
  }
  void range_add(int id, int l, int r, int ql, int qr, i64 v) {
    if (r < ql || l > qr) {
      return;
    }
    if (ql <= l && r <= qr) {
      apply_add(id, l, r, v);
      return;
    }
    push(id, l, r);
    int m = l + r >> 1;
    range_add(id << 1, l, m, ql, qr, v);
    range_add(id << 1 | 1, m + 1, r, ql, qr, v);
    pull(id);
  }
  void range_add(int ql, int qr, i64 v) {
    range_add(1, 0, n - 1, ql, qr, v);
  }
  i64 range_sum(int id, int l, int r, int ql, int qr) {
    if (r < ql || l > qr) {
      return 0;
    }
    if (ql <= l && r <= qr) {
      return t[id].sum;
    }
    push(id, l, r);
    int m = l + r >> 1;
    return range_sum(id << 1, l, m, ql, qr) + range_sum(id << 1 | 1, m + 1, r, ql, qr);
  }
  i64 range_sum(int ql, int qr) {
    return range_sum(1, 0, n - 1, ql, qr);
  }
};
