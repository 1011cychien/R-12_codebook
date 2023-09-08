// need to discretize
struct SuffixArray {
    int n;
    vector<int> sa, as, ha;
template <typename T>
    vector<int> sais(const T &s) {
        int n = s.size(), m = *max_element(s.begin(), s.end()) + 1;
        vector<int> pos(m + 1), f(n);
        for (auto ch : s) { pos[ch + 1]++; }
        for (int i = 0; i < m; i++) { pos[i + 1] += pos[i]; }
        for (int i = n - 2; i >= 0; i--) { f[i] = s[i] != s[i + 1] ? s[i] < s[i + 1] : f[i + 1]; }
        vector<int> x(m), sa(n);
        auto induce = [&](const vector<int> &ls) {
            fill(sa.begin(), sa.end(), -1);
            auto L = [&](int i) { if (i >= 0 && !f[i]) { sa[x[s[i]]++] = i; }};
            auto S = [&](int i) { if (i >= 0 && f[i]) { sa[--x[s[i]]] = i; }};
            for (int i = 0; i < m; i++) { x[i] = pos[i + 1]; }
            for (int i = int(ls.size()) - 1; i >= 0; i--) { S(ls[i]); }
            for (int i = 0; i < m; i++) { x[i] = pos[i]; }
            L(n - 1);
            for (int i = 0; i < n; i++) { L(sa[i] - 1); }
            for (int i = 0; i < m; i++) { x[i] = pos[i + 1]; }
            for (int i = n - 1; i >= 0; i--) { S(sa[i] - 1); }
        };
        auto ok = [&](int i) { return i == n || !f[i - 1] && f[i]; };
        auto same = [&](int i, int j) {
            do { if (s[i++] != s[j++]) { return false; }} while (!ok(i) && !ok(j));
            return ok(i) && ok(j);
        };
        vector<int> val(n), lms;
        for (int i = 1; i < n; i++) { if (ok(i)) { lms.push_back(i); }}
        induce(lms);
        if (!lms.empty()) {
            int p = -1, w = 0;
            for (auto v : sa) {
                if (v != 0 && ok(v)) {
                    if (p != -1 && same(p, v)) { w--; }
                    val[p = v] = w++;
                }
            }
            auto b = lms;
            for (auto &v : b) { v = val[v]; }
            b = sais(b);
            for (auto &v : b) { v = lms[v]; }
            induce(b);
        }
        return sa;
    }
template <typename T>
    SuffixArray(const T &s) : n(s.size()), sa(sais(s)), as(n), ha(n - 1) {
        for (int i = 0; i < n; i++) { as[sa[i]] = i; }
        for (int i = 0, j = 0; i < n; ++i) {
            if (as[i] == 0) {
                j = 0;
            } else {
                for (j -= j > 0; i + j < n && sa[as[i] - 1] + j < n && s[i + j] == s[sa[as[i] - 1] + j]; ) { ++j; }
                ha[as[i] - 1] = j;
            }
        }
    }
};