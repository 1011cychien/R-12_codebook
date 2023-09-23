// p : a[0] ~ a[d - 1]
// q : a[i] = \sum a[i - j]q[j]
template <typename T>
T linearRecurrence(vector<T> p, vector<T> q, i64 n) {
    int d = q.size() - 1;
    assert(int(p.size()) == d);
    p = p * q;
    p.resize(d);
    while (n > 0) {
        auto nq = q;
        for (int i = 1; i <= d; i += 2) {
            nq[i] *= -1;
        }
        auto np = p * nq;
        nq = q * nq;
        for (int i = 0; i < d; i++) {
            p[i] = np[i * 2 + n % 2];
        }
        for (int i = 0; i <= d; i++) {
            q[i] = nq[i * 2];
        }
        n /= 2;
    }
    return p[0] / q[0];
}