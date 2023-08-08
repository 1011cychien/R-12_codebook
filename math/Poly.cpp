template <int MOD>
struct Modint {
    static constexpr int P = MOD;
    int v;
    constexpr Modint() : v(0) {}
    constexpr Modint(i64 v_) : v(v_ % P) { if (v < 0) { v += P; }}
    constexpr friend Modint operator+(Modint a, Modint b) { return Modint((a.v + b.v) % P); }
    constexpr friend Modint operator-(Modint a, Modint b) { return Modint((a.v + P - b.v) % P); }
    constexpr friend Modint operator*(Modint a, Modint b) { return Modint(1LL * a.v * b.v % P); }
    constexpr Modint qpow(i64 p) {
        Modint res = 1, x = v;
        while (p > 0) {
            if (p & 1) { res = res * x; }
            x = x * x;
            p >>= 1;
        }
        return res;
    }
    constexpr Modint inv() { return qpow(P - 2); }
};
template<int P>
constexpr Modint<P> findPrimitiveRoot() {
    Modint<P> i = 2;
    int k = __builtin_ctz(P - 1);
    while (true) {
        if (i.qpow((P - 1) / 2).v != 1) { break; }
        i = i + 1;
    }
    return i.qpow(P - 1 >> k);
}
template<int P>
constexpr Modint<P> primitiveRoot = findPrimitiveRoot<P>();
vector<int> rev;
template<int P>
vector<Modint<P>> roots{0, 1};
template<int P>
void dft(vector<Modint<P>> &a) {
    int n = a.size();
    if (int(rev.size()) != n) {
        int k = __builtin_ctz(n) - 1;
        rev.resize(n);
        for (int i = 0; i < n; i++) { rev[i] = rev[i >> 1] >> 1 | (i & 1) << k; }
    }
    for (int i = 0; i < n; i++) { if (rev[i] < i) { swap(a[i], a[rev[i]]); }}
    if (roots<P>.size() < n) {
        int k = __builtin_ctz(roots<P>.size());
        roots<P>.resize(n);
        while ((1 << k) < n) {
            auto e = Modint<P>(primitiveRoot<P>).qpow(P - 1 >> k + 1);
            for (int i = 1 << k - 1; i < 1 << k; i++) {
                roots<P>[2 * i] = roots<P>[i];
                roots<P>[2 * i + 1] = roots<P>[i] * e;
            }
            k++;
        }
    }
    for (int k = 1; k < n; k *= 2) {
        for (int i = 0; i < n; i += 2 * k) {
            for (int j = 0; j < k; j++) {
                Modint<P> u = a[i + j];
                Modint<P> v = a[i + j + k] * roots<P>[k + j];
                a[i + j] = u + v;
                a[i + j + k] = u - v;
            }
        }
    }
}
template <int P>
void idft(vector<Modint<P>> &a) {
    int n = a.size();
    reverse(a.begin() + 1, a.end());
    dft(a);
    Modint<P> x = (1 - P) / n;
    for (int i = 0; i < n; i++) { a[i] = a[i] * x; }
}
template <int P>
struct Poly {
    using Mint = Modint<P>;
    vector<Mint> a;
    Poly() {}
    explicit Poly(const vector<Mint> &a) : a(a) {}
    explicit Poly(const initializer_list<Mint> &a) : a(a) {}
    explicit Poly(int n) : a(n) {}
template<class F>
    explicit Poly(int n, F f) : a(n) {
        for (int i = 0; i < n; i++) { a[i] = f(i); }
    }
    int size() const { return a.size(); }
    void resize(int n) { a.resize(n); }
    Mint operator[](int idx) const {
        if (idx < 0 || idx >= size()) { return 0; }
        return a[idx];
    }
    Mint& operator[](int idx) { return a[idx]; }
    Poly mulxk(int k) {
        auto b = a;
        b.insert(b.begin(), k, 0);
        return Poly(b);
    }
    Poly modxk(int k) {
        k = min(k, size());
        return Poly(vector<Mint>(a.begin(), a.begin() + k));
    }
    Poly divxk(int k) {
        if (size() <= k) { return Poly(); }
        return Poly(vector<Mint>(a.begin() + k, a.end()));
    }
    friend Poly operator+(const Poly &a, const Poly &b) {
        vector<Mint> res(max(a.size(), b.size()));
        for (int i = 0; i < int(res.size()); i++) { res[i] = a[i] + b[i]; }
        return Poly(res);
    }
    friend Poly operator-(const Poly &a, const Poly &b) {
        vector<Mint> res(max(a.size(), b.size()));
        for (int i = 0; i < int(res.size()); i++) { res[i] = a[i] - b[i]; }
        return Poly(res);
    }
    friend Poly operator*(Poly a, Poly b) {
        if (a.size() == 0 || b.size() == 0) { return Poly(); }
        int sz = 1, tot = a.size() + b.size() - 1;
        while (sz < tot) { sz *= 2; }
        a.resize(sz);
        b.resize(sz);
        dft(a.a);
        dft(b.a);
        for (int i = 0; i < sz; i++) { a[i] = a[i] * b[i]; }
        idft(a.a);
        a.resize(tot);
        return a;
    }
    friend Poly operator*(Mint a, Poly b) {
        for (int i = 0; i < int(b.size()); i++) { b[i] = b[i] * a; }
        return b;
    }
    friend Poly operator*(Poly a, Mint b) {
        for (int i = 0; i < int(a.size()); i++) { a[i] = a[i] * b; }
        return a;
    }
    Poly derivative()  {
        if (a.empty()) { return Poly(); }
        vector<Mint> res(size() - 1);
        for (int i = 0; i < size() - 1; ++i) { res[i] = (i + 1) * a[i + 1]; }
        return Poly(res);
    }
    Poly integral() {
        vector<Mint> res(size() + 1);
        for (int i = 0; i < size(); ++i) { res[i + 1] = a[i] * Mint(i + 1).inv(); }
        return Poly(res);
    }
    Poly inv(int m) {
        // a[0] != 0
        Poly x({a[0].inv()});
        int k = 1;
        while (k < m) {
            k *= 2;
            x = (x * (Poly({2}) - modxk(k) * x)).modxk(k);
        }
        return x.modxk(m);
    }
    Poly log(int m)  {
        return (derivative() * inv(m)).integral().modxk(m);
    }
    Poly exp(int m) {
        Poly x({1});
        int k = 1;
        while (k < m) {
            k *= 2;
            x = (x * (Poly({1}) - x.log(k) + modxk(k))).modxk(k);
        }
        return x.modxk(m);
    }
    Poly pow(i64 k, int m) {
        if (k == 0) {
            vector<Mint> x(m);
            x[0] = 1;
            return Poly(x);
        }
        int i = 0;
        while (i < size() && a[i].v == 0) { i++; }
        if (i == size() || __int128(i) * k >= m) { return Poly(vector<Mint>(m)); }
        Mint v = a[i];
        auto f = divxk(i) * v.inv();
        return (f.log(m - i * k) * k).exp(m - i * k).mulxk(i * k) * v.qpow(k);
    }
    Poly sqrt(int m) {
        // a[0] == 1, otherwise quadratic residue?
        Poly x({1});
        int k = 1;
        while (k < m) {
            k *= 2;
            x = (x + (modxk(k) * x.inv(k)).modxk(k)) * ((Mint::P + 1) / 2);
        }
        return x.modxk(m);
    }
    Poly mulT(Poly b) const {
        if (b.size() == 0) { return Poly(); }
        int n = b.size();
        reverse(b.a.begin(), b.a.end());
        return ((*this) * b).divxk(n - 1);
    }
    vector<Mint> evaluation(vector<Mint> x) {
        if (size() == 0) { return vector<Mint>(x.size(), 0); }
        const int n = max(int(x.size()), size());
        vector<Poly> q(4 * n);
        vector<Mint> ans(x.size());
        x.resize(n);
        auto build = [&](auto build, int id, int l, int r) -> void {
            if (r - l == 1) {
                q[id] = Poly({1, -x[l].v});
            } else {
                int m = (l + r) / 2;
                build(build, 2 * id, l, m);
                build(build, 2 * id + 1, m, r);
                q[id] = q[2 * id] * q[2 * id + 1];
            }
        };
        build(build, 1, 0, n);
        auto work = [&](auto work, int id, int l, int r, const Poly &num) -> void {
            if (r - l == 1) {
                if (l < int(ans.size())) { ans[l] = num[0]; }
            } else {
                int m = (l + r) / 2;
                work(work, 2 * id, l, m, num.mulT(q[2 * id + 1]).modxk(m - l));
                work(work, 2 * id + 1, m, r, num.mulT(q[2 * id]).modxk(r - m));
            }
        };
        work(work, 1, 0, n, mulT(q[1].inv(n)));
        return ans;
    }
};
template <int P>
Poly<P> interpolate(vector<Modint<P>> x, vector<Modint<P>> y) {
    // f(xi) = yi
    int n = x.size();
    vector<Poly<P>> p(4 * n), q(4 * n);
    auto dfs1 = [&](auto dfs1, int id, int l, int r) -> void {
        if (l == r) {
            p[id] = Poly<P>({-x[l].v, 1});
            return;
        }
        int m = l + r >> 1;
        dfs1(dfs1, id << 1, l, m);
        dfs1(dfs1, id << 1 | 1, m + 1, r);
        p[id] = p[id << 1] * p[id << 1 | 1];
    };
    dfs1(dfs1, 1, 0, n - 1);
    Poly<P> f = Poly<P>(p[1].derivative().evaluation(x));
    auto dfs2 = [&](auto dfs2, int id, int l, int r) -> void {
        if (l == r) {
            q[id] = Poly<P>({y[l] * f[l].inv()});
            return;
        }
        int m = l + r >> 1;
        dfs2(dfs2, id << 1, l, m);
        dfs2(dfs2, id << 1 | 1, m + 1, r);
        q[id] = q[id << 1] * p[id << 1 | 1] + q[id << 1 | 1] * p[id << 1];
    };
    dfs2(dfs2, 1, 0, n - 1);
    return q[1];
}