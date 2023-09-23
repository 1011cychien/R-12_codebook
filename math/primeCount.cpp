// __attribute__((target("avx2"), optimize("O3", "unroll-loops")))
i64 primeCount(const i64 n) {
    if (n <= 1) { return 0; }
    if (n == 2) { return 1; }
    const int v = sqrtl(n);
    int s = (v + 1) / 2;
    vector<int> smalls(s), roughs(s), skip(v + 1);
    vector<i64> larges(s);
    iota(smalls.begin(), smalls.end(), 0);
    for (int i = 0; i < s; i++) {
        roughs[i] = 2 * i + 1;
        larges[i] = (n / roughs[i] - 1) / 2;
    }
    const auto half = [](int n) -> int { return (n - 1) >> 1; };
    int pc = 0;
    for (int p = 3; p <= v; p += 2) {
        if (skip[p]) { continue; }
        int q = p * p;
        if (1LL * q * q > n) { break; }
        skip[p] = true;
        for (int i = q; i <= v; i += 2 * p) skip[i] = true;
        int ns = 0;
        for (int k = 0; k < s; k++) {
            int i = roughs[k];
            if (skip[i]) { continue; }
            i64 d = 1LL * i * p;
            larges[ns] = larges[k] - (d <= v ? larges[smalls[d / 2] - pc] : smalls[half(n / d)]) + pc;
            roughs[ns++] = i;
        }
        s = ns;
        for (int i = half(v), j = v / p - 1 | 1; j >= p; j -= 2) {
            int c = smalls[j / 2] - pc;
            for (int e = j * p / 2; i >= e; i--) { smalls[i] -= c; }
        }
        pc++;
    }
    larges[0] += 1LL * (s + 2 * (pc - 1)) * (s - 1) / 2;
    for (int k = 1; k < s; k++) { larges[0] -= larges[k]; }
    for (int l = 1; l < s; l++) {
        i64 q = roughs[l];
        i64 M = n / q;
        int e = smalls[half(M / q)] - pc;
        if (e <= l) { break; }
        i64 t = 0;
        for (int k = l + 1; k <= e; k++) { t += smalls[half(M / roughs[k])]; }
        larges[0] += t - 1LL * (e - l) * (pc + l - 1);
    }
    return larges[0] + 1;
}