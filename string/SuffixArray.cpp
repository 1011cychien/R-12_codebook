struct SuffixArray {
    int n;
    vector<int> sa, as, ha;
    vector<vector<int>> rmq;
template <typename T>
    SuffixArray(const T &s) : n(s.size()), sa(n), as(n), ha(n - 1) {
        n = s.size();
        iota(sa.begin(), sa.end(), 0);
        sort(sa.begin(), sa.end(), [&](int a, int b) { return s[a] < s[b]; });
        as[sa[0]] = 0;
        for (int i = 1; i < n; ++i) {
            as[sa[i]] = as[sa[i - 1]] + (s[sa[i]] != s[sa[i - 1]]);
        }
        int k = 1;
        vector<int> tmp, cnt(n);
        tmp.reserve(n);
        while (as[sa[n - 1]] < n - 1) {
            tmp.clear();
            for (int i = 0; i < k; ++i) { tmp.push_back(n - k + i); }
            for (auto i : sa) { if (i >= k) { tmp.push_back(i - k); } }
            fill(cnt.begin(), cnt.end(), 0);
            for (int i = 0; i < n; ++i) { ++cnt[as[i]]; }
            for (int i = 1; i < n; ++i) { cnt[i] += cnt[i - 1]; }
            for (int i = n - 1; i >= 0; --i) { sa[--cnt[as[tmp[i]]]] = tmp[i]; }
            swap(as, tmp);
            as[sa[0]] = 0;
            for (int i = 1; i < n; ++i) {
                as[sa[i]] = as[sa[i - 1]] + (tmp[sa[i - 1]] < tmp[sa[i]] || sa[i - 1] + k == n || tmp[sa[i - 1] + k] < tmp[sa[i] + k]);
            }
            k *= 2;
        }
        for (int i = 0, j = 0; i < n; ++i) {
            if (as[i] == 0) {
                j = 0;
            } else {
                for (j -= j > 0; i + j < n && sa[as[i] - 1] + j < n && s[i + j] == s[sa[as[i] - 1] + j]; ) { ++j; }
                ha[as[i] - 1] = j;
            }
        }
        if (n > 1) {
            const int lg = __lg(n - 1) + 1;
            rmq.assign(lg + 1, vector<int>(n - 1));
            rmq[0] = ha;
            for (int i = 1; i <= lg; i++) {
                for (int j = 0; j + (1 << i) <= n; j++) {
                    rmq[i][j] = min(rmq[i - 1][j], rmq[i - 1][j + (1 << i - 1)]);
                }
            }
        }
    }
    int lcp(int x, int y) {
        if (x == y) { return n - x; }
        x = as[x], y = as[y];
        if (x > y) { swap(x, y); }
        int k = __lg(y - x);
        return min(rmq[k][x], rmq[k][y - (1 << k)]);
    }
};