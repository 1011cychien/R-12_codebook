// signed, may not be integer
template <typename T>
dbl area(vector<P<T>> a) {
  T res = 0;
  int n = a.size();
  for (int i = 0; i < n; i++) {
    res += cross(a[i], a[(i + 1) % n]);
  }
  return res / 2.0;
}

template <typename T>
vector<P<T>> convexHull(vector<P<T>> a) {
  int n = a.size();
  if (n <= 1) {
    return a;
  }
  sort(a.begin(), a.end());
  vector<P<T>> res(2 * n);
  int j = 0;
  for (int i = 0; i < n; res[j++] = a[i++]) {
    while (j >= 2 && side(res[j - 2], res[j - 1], a[i]) <= 0) { j--; }
  }
  for (int i = n - 2, k = j; i >= 0; res[j++] = a[i--]) {
    while (j > k && side(res[j - 2], res[j - 1], a[i]) <= 0) { j--; }
  }
  res.resize(j - 1);
  return res;
}
