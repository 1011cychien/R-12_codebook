// \min i : [0, n) (a * i + b) % m
// ok in 1e9
int minModLinear(int n, int m, int a, int b, int cnt = 1, int p = 1, int q = 1) {
  if (a == 0) { return b; }
  if (cnt % 2 == 1) {
    if (b >= a) {
      int t = (m - b + a - 1) / a;
      int c = (t - 1) * p + q;
      if (n <= c) { return b; }
      n -= c;
      b += a * t - m;
    }
    b = a - 1 - b;
  } else {
    if (b < m - a) {
      int t = (m - b - 1) / a;
      int c = t * p;
      if (n <= c) { return (n - 1) / p * a + b; }
      n -= c;
      b += a * t;
    }
    b = m - 1 - b;
  }
  cnt++;
  int d = m / a;
  int c = minModLinear(n, a, m % a, b, cnt, (d - 1) * p + q, d * p + q);
  return cnt % 2 == 1 ? m - 1 - c : a - 1 - c;
}
