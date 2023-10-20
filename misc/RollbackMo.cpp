for (int l = 0, r = -1; auto [ql, qr, i] : qs) {
  if (ql / B == qr / B) {
    for (int j = ql; j <= qr; j++) {
      cntSmall[a[j]]++;
      ans[i] = max(ans[i], 1LL * b[a[j]] * cntSmall[a[j]]);
    }
    for (int j = ql; j <= qr; j++) {
      cntSmall[a[j]]--;
    }
    continue;
  }
  if (int block = ql / B; block != lst) {
    int x = min((block + 1) * B, n);
    while (r + 1 < x) { add(++r); }
    while (r >= x) { del(r--); }
    while (l < x) { del(l++); }
    mx = 0;
    lst = block;
  }
  while (r < qr) { add(++r); }
  i64 tmpMx = mx;
  int tmpL = l;
  while (l > ql) { add(--l); }
  ans[i] = mx;
  mx = tmpMx;
  while (l < tmpL) { del(l++); }
}
