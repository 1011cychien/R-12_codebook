vector<int> euler, vis(V);
auto dfs = [&](auto dfs, int u) -> void {
    while (!adj[u].empty()) {
        while (!adj[u].empty() && del[adj[u].back()[1]]) {
            adj[u].pop_back();
        }
        if (!adj[u].empty()) {
            auto [v, i] = adj[u].back();
            del[i] = true;
            dfs(dfs, v);
        }
    }
    euler.push_back(u);
};            
dfs(dfs, 0);
reverse(euler.begin(), euler.end());