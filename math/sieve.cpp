vector<int> minp(N + 1), primes, mobius(N + 1);
mobius[1] = 1;
for (int i = 2; i <= N; i++) {
  if (!minp[i]) {
    primes.push_back(i);
    minp[i] = i;
    mobius[i] = -1;
  }
  for (int p : primes) {
    if (p > N / i) {
      break;
    }
    minp[p * i] = p;
    mobius[p * i] = -mobius[i];
    if (i % p == 0) {
      mobius[p * i] = 0;
      break;
    }
  }
}
