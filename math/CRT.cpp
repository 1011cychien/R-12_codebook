// returns (rem, mod), n = 0 return (0, 1), no solution return (0, 0)
pair<i64, i64> crt(vector<i64> r, vector<i64> m) {
    int n = r.size();
    for (int i = 0; i < n; i++) {
        r[i] %= m[i];
        if (r[i] < 0) { r[i] += m[i]; }
    }
    i64 r0 = 0, m0 = 1;
    for (int i = 0; i < n; i++) {
        i64 r1 = r[i], m1 = m[i];
        if (m0 < m1) { swap(r0, r1), swap(m0, m1); }
        if (m0 % m1 == 0) {
            if (r0 % m1 != r1) { return {0, 0}; }
            continue;
        }
        auto [g, a, b] = extgcd(m0, m1);
        i64 u1 = m1 / g;
        if ((r1 - r0) % g != 0) { return {0, 0}; }
        i64 x = (r1 - r0) / g % u1 * a % u1;
        r0 += x * m0;
        m0 *= u1;
        if (r0 < 0) { r0 += m0; }
    }
    return {r0, m0};
}