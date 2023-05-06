auto cmp = [&](int i, int j) { return tin[i] < tin[j]; };
sort(verts.begin(), verts.end(), cmp);
for (int i = int(verts.size()) - 1; i > 0; i--) {
    verts.push_back(lca(verts[i], verts[i - 1]));
}
sort(verts.begin(), verts.end(), cmp);
verts.erase(unique(verts.begin(), verts.end()), verts.end());
vector<int> stk;
for (auto u : verts) {
    G[u].clear();
    while (!stk.empty() && tin[u] > tout[stk.back()]) {
        stk.pop_back();
    }
    if (!stk.empty()) {
        G[stk.back()].push_back(u);
    }
    stk.push_back(u);
}