constexpr int P = 998244353, RT = 3;
// qpow(int x, i64 p) 
vector<int> rev;
vector<int> roots{0, 1};
void dft(vector<int> &a) {
    int n = a.size();
    if (int(rev.size()) != n) {
        int k = __builtin_ctz(n) - 1;
        rev.resize(n);
        for (int i = 0; i < n; i++) {
            rev[i] = rev[i >> 1] >> 1 | (i & 1) << k;
        }
    }
    for (int i = 0; i < n; i++) {
        if (rev[i] < i) {
            swap(a[i], a[rev[i]]);
        }
    }
    if (int(roots.size()) < n) {
        int k = __builtin_ctz(roots.size());
        roots.resize(n);
        while ((1 << k) < n) {
            int e = qpow(RT, P - 1 >> k + 1);
            for (int i = 1 << k - 1; i < 1 << k; i++) {
                roots[2 * i] = roots[i];
                roots[2 * i + 1] = 1LL * roots[i] * e % P;
            }
            k++;
        }
    }
    for (int k = 1; k < n; k *= 2) {
        for (int i = 0; i < n; i += 2 * k) {
            for (int j = 0; j < k; j++) {
                int u = a[i + j], v = 1LL * a[i + j + k] * roots[k + j] % P;
                a[i + j] = (u + v) % P;
                a[i + j + k] = (u + P - v) % P;
            }
        }
    }
}
void idft(vector<int> &a) {
    int n = a.size();
    reverse(a.begin() + 1, a.end());
    dft(a);
    int x = P + (1 - P) / n;
    for (int i = 0; i < n; i++) {
        a[i] = 1LL * a[i] * x % P;
    }
}

struct Poly {
    vector<int> a;
    Poly() {}
    explicit Poly(const vector<int> &a) : a(a) {}
    explicit Poly(const initializer_list<int> &a) : a(a) {}
    explicit Poly(int n) : a(n) {}
template<class F>
    explicit Poly(int n, F f) : a(n) {
        for (int i = 0; i < n; i++) {
            a[i] = f(i);
        }
    }
    int size() const {
        return a.size();
    }
    void resize(int n) {
        a.resize(n);
    }
    int operator[](int idx) const {
        if (idx < 0 || idx >= size()) {
            return 0;
        }
        return a[idx];
    }
    int& operator[](int idx) {
        return a[idx];
    }
    Poly mulxk(int k) const {
        auto b = a;
        b.insert(b.begin(), k, 0);
        return Poly(b);
    }
    Poly modxk(int k) const {
        k = min(k, size());
        return Poly(vector<int>(a.begin(), a.begin() + k));
    }
    Poly divxk(int k) const {
        if (size() <= k) {
            return Poly();
        }
        return Poly(vector<int>(a.begin() + k, a.end()));
    }
    friend Poly operator+(const Poly &a, const Poly &b) {
        vector<int> res(max(a.size(), b.size()));
        for (int i = 0; i < int(res.size()); i++) {
            res[i] = (a[i] + b[i]) % P;
        }
        return Poly(res);
    }
    friend Poly operator-(const Poly &a, const Poly &b) {
        vector<int> res(max(a.size(), b.size()));
        for (int i = 0; i < int(res.size()); i++) {
            res[i] = (a[i] + P - b[i]) % P;
        }
        return Poly(res);
    }
    friend Poly operator*(Poly a, Poly b) {
        if (a.size() == 0 || b.size() == 0) {
            return Poly();
        }
        int sz = 1, tot = a.size() + b.size() - 1;
        while (sz < tot) { sz *= 2; }
        a.resize(sz);
        b.resize(sz);
        dft(a.a);
        dft(b.a);
        for (int i = 0; i < sz; i++) {
            a[i] = 1LL * a[i] * b[i] % P;
        }
        idft(a.a);
        a.resize(tot);
        return a;
    }
    friend Poly operator*(i64 a, Poly b) {
        for (int i = 0; i < int(b.size()); i++) {
            b[i] = a % P * b[i] % P;
        }
        return b;
    }
    friend Poly operator*(Poly a, i64 b) {
        for (int i = 0; i < int(a.size()); i++) {
            a[i] = b % P * a[i] % P;
        }
        return a;
    }
    Poly& operator+=(Poly b) {
        return (*this) = (*this) + b;
    }
    Poly& operator-=(Poly b) {
        return (*this) = (*this) - b;
    }
    Poly& operator*=(Poly b) {
        return (*this) = (*this) * b;
    }
    Poly derivative() const {
        if (a.empty()) { return Poly(); }
        vector<int> res(size() - 1);
        for (int i = 0; i < size() - 1; ++i) {
            res[i] = 1LL * (i + 1) * a[i + 1] % P;
        }
        return Poly(res);
    }
    Poly integral() const {
        vector<int> res(size() + 1);
        for (int i = 0; i < size(); ++i) {
            res[i + 1] = 1LL * a[i] * qpow(i + 1, P - 2) % P;
        }
        return Poly(res);
    }
    Poly inv(int m) const {
        Poly x({qpow(a[0], P - 2)});
        int k = 1;
        while (k < m) {
            k *= 2;
            x = (x * (Poly({2}) - modxk(k) * x)).modxk(k);
        }
        return x.modxk(m);
    }
    Poly log(int m) const {
        return (derivative() * inv(m)).integral().modxk(m);
    }
    Poly exp(int m) const {
        Poly x({1});
        int k = 1;
        while (k < m) {
            k *= 2;
            x = (x * (Poly({1}) - x.log(k) + modxk(k))).modxk(k);
        }
        return x.modxk(m);
    }
    Poly pow(i64 k, int m) const {
        if (k == 0) {
            return Poly(m, [&](int i) { return i == 0; });
        }
        int i = 0;
        while (i < size() && a[i] == 0) {
            i++;
        }
        if (i == size() || __int128(i) * k >= m) {
            return Poly(m);
        }
        int v = a[i];
        auto f = divxk(i) * qpow(v, P - 2);
        return (f.log(m - i * k) * k).exp(m - i * k).mulxk(i * k) * qpow(v, k);
    }
    Poly sqrt(int m) const {
        // a[0] = 1
        Poly x({1});
        int k = 1;
        while (k < m) {
            k *= 2;
            x = (x + (modxk(k) * x.inv(k)).modxk(k)) * ((P + 1) / 2);
        }
        return x.modxk(m);
    }
    Poly mulT(Poly b) const {
        if (b.size() == 0) {
            return Poly();
        }
        int n = b.size();
        reverse(b.a.begin(), b.a.end());
        return ((*this) * b).divxk(n - 1);
    }
    vector<int> evaluation(vector<int> x) const {
        if (size() == 0) {
            return vector<int>(x.size());
        }
        const int n = max(int(x.size()), size());
        vector<Poly> q(4 * n);
        vector<int> ans(x.size());
        x.resize(n);
        auto build = [&](auto build, int p, int l, int r) -> void {
            if (r - l == 1) {
                q[p] = Poly({1, (P - x[l]) % P});
            } else {
                int m = (l + r) / 2;
                build(build, 2 * p, l, m);
                build(build, 2 * p + 1, m, r);
                q[p] = q[2 * p] * q[2 * p + 1];
            }
        };
        build(build, 1, 0, n);
        auto work = [&](auto work, int p, int l, int r, const Poly &num) -> void {
            if (r - l == 1) {
                if (l < int(ans.size())) {
                    ans[l] = num[0];
                }
            } else {
                int m = (l + r) / 2;
                work(work, 2 * p, l, m, num.mulT(q[2 * p + 1]).modxk(m - l));
                work(work, 2 * p + 1, m, r, num.mulT(q[2 * p]).modxk(r - m));
            }
        };
        work(work, 1, 0, n, mulT(q[1].inv(n)));
        return ans;
    }
};
vector<int> interpolate(vector<int> x, vector<int> y) {
    // f(xi) = yi
    int n = x.size();
    vector<Poly> p(4 * n), q(4 * n);
    auto dfs1 = [&](auto dfs1, int id, int l, int r) -> void {
        if (l == r) {
            p[id] = Poly({(P - x[l]) % P, 1});
            return;
        }
        int m = l + r >> 1;
        dfs1(dfs1, id << 1, l, m);
        dfs1(dfs1, id << 1 | 1, m + 1, r);
        p[id] = p[id << 1] * p[id << 1 | 1];
    };
    dfs1(dfs1, 1, 0, n - 1);
    Poly f = Poly(p[1].derivative().evaluation(x));
    auto dfs2 = [&](auto dfs2, int id, int l, int r) -> void {
        if (l == r) {
            q[id] = Poly({int(1LL * y[l] * qpow(f[l], P - 2) % P)});
            return;
        }
        int m = l + r >> 1;
        dfs2(dfs2, id << 1, l, m);
        dfs2(dfs2, id << 1 | 1, m + 1, r);
        q[id] = q[id << 1] * p[id << 1 | 1] + q[id << 1 | 1] * p[id << 1];
    };
    dfs2(dfs2, 1, 0, n - 1);
    return q[1].a;
}