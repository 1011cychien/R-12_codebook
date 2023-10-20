int n, T;
cin >> n >> T;
vector<int> cnt(T + 1);
for (int i = 0; i < n; i++) {
  int a;
  cin >> a;
  cnt[a]++;
}
vector<Mint> inv(T + 1);
for (int i = 1; i <= T; i++) {
  inv[i] = i == 1 ? 1 : -P / i * inv[P % i];
}
FPS f(T + 1);
for (int i = 1; i <= T; i++) {
  for (int j = 1; j * i <= T; j++) {
    f[i * j] = f[i * j] + (j % 2 == 1 ? 1 : -1) * cnt[i] * inv[j];
  }
}
f = f.exp(T + 1);
