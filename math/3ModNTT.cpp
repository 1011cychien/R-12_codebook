constexpr int P0 = 998244353, P1 = 1004535809, P2 = 469762049;
constexpr i64 P01 = 1LL * P0 * P1;
constexpr int inv0 = Modint<P1>(P0).inv().v;
constexpr int inv01 = Modint<P2>(P01).inv().v;
for (int i = 0; i < int(c.size()); i++) {
  i64 x = 1LL * (c1[i] - c0[i] + P1) % P1 * inv0 % P1 * P0 + c0[i];
  c[i] = ((c2[i] - x % P2 + P2) % P2 * inv01 % P2 * (P01 % P) % P + x) % P;
}
