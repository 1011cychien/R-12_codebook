// bipartite
vector<int> ans(m, -1);
vector has(a + b, vector<pair<int, int>>(col, {-1, -1}));
for (int i = 0; i < m; i++) {
    auto [u, v] = e[i];
    vector<int> c;
    for (auto x : {u, v}) {
        c.push_back(0);
        while (has[x][c.back()].first != -1) {
            c.back()++;
        }
    }
    if (c[0] != c[1]) {
        auto dfs = [&](auto self, int u, int x) -> void {
            auto [v, i] = has[u][c[x]];
            if (v != -1) {
                if (has[v][c[x ^ 1]].first != -1) {
                    self(self, v, x ^ 1);
                } else {
                    has[v][c[x]] = {-1, -1};
                }
                has[u][c[x ^ 1]] = {v, i}, has[v][c[x ^ 1]] = {u, i};
                ans[i] = c[x ^ 1];
            }
        };
        dfs(dfs, v, 0);
    }
    has[u][c[0]] = {v, i};
    has[v][c[0]] = {u, i};
    ans[i] = c[0];
}
// general
auto vizing(int n, const vector<pair<int, int>> &e) {
    vector<int> deg(n);
    for (auto [u, v] : e) {
        deg[u]++, deg[v]++;
    }
    int col = *max_element(deg.begin(), deg.end()) + 1;
    vector<int> free(n);
    vector ans(n, vector<int>(n, -1));
    vector at(n, vector<int>(col, -1));
    auto update = [&](int u) {
        free[u] = 0;
        while (at[u][free[u]] != -1) {
            free[u]++;
        }
    };
    auto color = [&](int u, int v, int c1) {
        int c2 = ans[u][v];
        ans[u][v] = ans[v][u] = c1;
        at[u][c1] = v, at[v][c1] = u;
        if (c2 != -1) {
            at[u][c2] = at[v][c2] = -1;
            free[u] = free[v] = c2;
        } else {
            update(u), update(v);
        }
        return c2;
    };
    auto flip = [&](int u, int c1, int c2) {
        int v = at[u][c1];
        swap(at[u][c1], at[u][c2]);
        if (v != -1) {
            ans[u][v] = ans[v][u] = c2;
        }
        if (at[u][c1] == -1) {
            free[u] = c1;
        }
        if (at[u][c2] == -1) {
            free[u] = c2;
        }
        return v;
    };
    for (int i = 0; i < int(e.size()); i++) {
        auto [u, v1] = e[i];
        int v2 = v1, c1 = free[u], c2 = c1, d;
        vector<pair<int, int>> fan;
        vector<int> vis(col);
        while (ans[u][v1] == -1) {
            fan.emplace_back(v2, d = free[v2]);
            if (at[v2][c2] == -1) {
                for (int j = int(fan.size()) - 1; j >= 0; j--) {
                    c2 = color(u, fan[j].first, c2);
                }
            } else if (at[u][d] == -1) {
                for (int j = int(fan.size()) - 1; j >= 0; j--) {
                    color(u, fan[j].first, fan[j].second);
                }
            } else if (vis[d] == 1) {
                break;
            } else {
                vis[d] = 1, v2 = at[u][d];
            }
        }
        if (ans[u][v1] == -1) {
            while (v2 != -1) {
                v2= flip(v2, c2, d);
                swap(c2, d);
            }
            if (at[u][c1] != -1) {
                int j = int(fan.size()) - 2;
                while (j >= 0 && fan[j].second != c2) {
                    j--;
                }
                while (j >= 0) {
                    color(u, fan[j].first, fan[j].second);
                    j--;
                }
            } else {
                i--;
            }
        }
    }
    return pair(col, ans);
}