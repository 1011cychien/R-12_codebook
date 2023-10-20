// \sum {i = 0} {n} floor((a * i + b) / c)
i64 floorSum(i64 a, i64 b, i64 c, i64 n) {
  if (n < 0) { return 0; }
  if (n == 0) { return b / c; }
  if (a == 0) { return b / c * (n + 1); }
  i64 res = 0;
  if (a >= c) { res += a / c * n * (n + 1) / 2, a %= c; }
  if (b >= c) { res += b / c * (n + 1), b %= c; }
  i64 m = (a * n + b) / c;
  return res + n * m - (m == 0 ? 0 : floorSum(c, c - b - 1, a, m - 1));
}
