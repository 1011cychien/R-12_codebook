vector<vector<array<int, 2>>> g(m + 2);
for (int i = 0; i <= m; i++) {
    if (i == m || used[i]) {
        auto tmpDsu = curDsu;
        for (int j = 0; auto [w, u, v] : edges) {
            if (i != j && used[j]) {
                tmpDsu.join(u, v);
            }
            j++;
        }
        for (int j = 0; auto [w, u, v] : edges) {
            if (!used[j] && !tmpDsu.same(u, v)) {
                g[i].push_back({j, w}); // i == m, S = m
            }
            j++;
        }
    }
}
for (int i = 0; i <= m; i++) {
    if (i == m || used[i]) {
        auto tmpDeg = curDeg;
        for (int j = 0; auto [w, u, v] : edges) {
            if (i != j && used[j] && u < k) {
                tmpDeg[u]++;
            }
            j++;
        }
        for (int j = 0; auto [w, u, v] : edges) {
            if (!used[j] && (u >= k || tmpDeg[u] < deg[u])) {
                g[j].push_back({i == m ? m + 1 : i, i == m ? 0 : -edges[i][0]});
            }
            j++;
        }
    }
}
vector<int> q{m}, inq(m + 2), dis(m + 2, 1e9), pre(m + 2, -1);
// spfa
int u = m + 1;
while (true) {
    u = pre[u];
    if (u == m) {
        break;
    }
    used[u] ^= 1;
}
res += dis[m + 1];