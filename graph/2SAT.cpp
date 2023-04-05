void add(int u, bool x) {
    g[2 * u + !x].push_back(2 * u + x);
}
void add_clause(int u, bool x, int v, bool y) {
    g[2 * u + !x].push_back(2 * v + y);
    g[2 * v + !y].push_back(2 * u + x);
}
// build scc
vector<int> ans(n);
for (int i = 0; i < n; i++) {
    if (scc.id[2 * i] == scc.id[2 * i + 1]) {
        break;
    }
    ans[i] = scc.id[2 * i] > scc.id[2 * i + 1];
}