pair<int, vector<int>> maxClique(const vector<bitset<N>> adj) {
    int n = adj.size();
    int mx = 0;
    vector<int> ans, cur;
    auto rec = [&](auto rec, bitset<N> s) -> void {
        int sz = s.count();
        if (int(cur.size()) > mx) { mx = cur.size(), ans = cur; }
        if (int(cur.size()) + sz <= mx) { return; }
        int e1 = -1, e2 = -1;
        vector<int> d(n);
        for (int i = 0; i < n; i++) {
            if (s[i]) {
                d[i] = (adj[i] & s).count();
                if (e1 == -1 || d[i] > d[e1]) { e1 = i; }
                if (e2 == -1 || d[i] < d[e2]) { e2 = i; }
            }
        }
        if (d[e1] >= sz - 2) {
            cur.push_back(e1);
            auto s1 = adj[e1] & s;
            rec(rec, s1);
            cur.pop_back();
            return;
        }
        cur.push_back(e2);
        auto s2 = adj[e2] & s;
        rec(rec, s2);
        cur.pop_back();
        s.reset(e2);
        rec(rec, s);
    };
    bitset<N> all;
    for (int i = 0; i < n; i++) {
        all.set(i);
    }
    rec(rec, all);
    return pair(mx, ans);
}