i64 mul(i64 a, i64 b, i64 mod) {}
i64 qpow(i64 x, i64 p, i64 mod) {}
bool isPrime(i64 n) {
  if (n == 1) { return false; }
  int r = __builtin_ctzll(n - 1);
  i64 d = n - 1 >> r;
  auto checkComposite = [&](i64 p) {
    i64 x = qpow(p, d, n);
    if (x == 1 || x == n - 1) { return false; }
    for (int i = 1; i < r; i++) {
      x = mul(x, x, n);
      if (x == n - 1) { return false; }
    }
    return true;
  };
  for (auto p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
    if (n == p) {
      return true;
    } else if (checkComposite(p)) {
      return false;
    }
  }
  return true;
}
vector<i64> pollardRho(i64 n) {
  vector<i64> res;
  auto work = [&](auto work, i64 n) {
    if (n <= 10000) {
      for (int i = 2; i * i <= n; i++) {
        while (n % i == 0) {
          res.push_back(i);
          n /= i;
        }
      }
      if (n > 1) { res.push_back(n); }
      return;
    } else if (isPrime(n)) {
      res.push_back(n);
      return;
    }
    i64 x0 = 2;
    auto f = [&](i64 x) { return (mul(x, x, n) + 1) % n; };
    while (true) {
      i64 x = x0, y = x0, d = 1, power = 1, lam = 0, v = 1;
      while (d == 1) {
        y = f(y);
        ++lam;
        v = mul(v, abs(x - y), n);
        if (lam % 127 == 0) {
          d = gcd(v, n);
          v = 1;
        }
        if (power == lam) {
          x = y;
          power *= 2;
          lam = 0;
          d = gcd(v, n);
          v = 1;
        }
      }
      if (d != n) {
        work(work, d);
        work(work, n / d);
        return;
      }
      ++x0;
    }
  };
  work(work, n);
  sort(res.begin(), res.end());
  return res;
}
