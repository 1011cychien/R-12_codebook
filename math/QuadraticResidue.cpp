// rng
int jacobi(int a, int m) {
    int s = 1;
    while (m > 1) {
        a %= m;
        if (a == 0) { return 0; }
        int r = __builtin_ctz(a);
        if (r % 2 == 1 && (m + 2 & 4) != 0) { s = -s; }
        a >>= r;
        if ((a & m & 2) != 0) { s = -s; }
        swap(a, m);
    }
    return s;
}
int quadraticResidue(int a, int p) {
    if (p == 2) { return a % 2; }
    int j = jacobi(a, p);
    if (j == 0 || j == -1) { return j; }
    int b, d;
    while (true) {
        b = rng() % p;
        d = (1LL * b * b + p - a) % p;
        if (jacobi(d, p) == -1) { break; }
    }
    int f0 = b, f1 = 1, g0 = 1, g1 = 0, tmp;
    for (int e = p + 1 >> 1; e > 0; e >>= 1) {
        if (e % 2 == 1) {
            tmp = (1LL * g0 * f0 + 1LL * d * g1 % p * f1 % p) % p;
            g1 = (1LL * g0 * f1 + 1LL * g1 * f0) % p;
            g0 = tmp;
        }
        tmp = (1LL * f0 * f0 + 1LL * d * f1 % p * f1 % p) % p;
        f1 = 2LL * f0 * f1 % p;
        f0 = tmp;
    }
    return g0;
}