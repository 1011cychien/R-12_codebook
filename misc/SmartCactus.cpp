auto work = [&](const vector<int> cycle) {
  // merge cycle info to u?
  int len = cycle.size(), u = cycle[0];
};
auto dfs = [&](auto dfs, int u, int p) {
  par[u] = p;
  vis[u] = 1;
  for (auto v : adj[u]) {
    if (v == p) { continue; }
    if (vis[v] == 0) {
      dfs(dfs, v, u);
      if (!cyc[v]) { // merge dp }
    } else if (vis[v] == 1) {
      for (int w = u; w != v; w = par[w]) {
        cyc[w] = 1;
      }
    } else {
      vector<int> cycle = {u};
      for (int w = v; w != u; w = par[w]) { cycle.push_back(w); }
      work(cycle);
    }
  }

  vis[u] = 2;
};
