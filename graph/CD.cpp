vector<int> sz(n), vis(n);
auto build = [&](auto build, int u, int p) -> void {
  sz[u] = 1;
  for (auto v : g[u]) {
    if (v != p && !vis[v]) {
      build(build, v, u);
      sz[u] += sz[v];
    }
  }
};
auto find = [&](auto find, int u, int p, int tot) -> int {
  for (auto v : g[u]) {
    if (v != p && !vis[v] && 2 * sz[v] > tot) {
      return find(find, v, u, tot);
    }
  }
  return u;
};

auto dfs = [&](auto dfs, int cen) -> void {
  build(build, cen, -1);
  cen = find(find, cen, -1, sz[cen]);
  vis[cen] = 1;
  build(build, cen, -1);

  for (auto v : g[cen]) {
    if (!vis[v]) {
      dfs(dfs, v);
    }
  }
};
dfs(dfs, 0);
