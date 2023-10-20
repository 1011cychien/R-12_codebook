sort(ord.begin(), ord.end(), [&](auto i, auto j) { return pair(deg[i], i) > pair(deg[j], j); });
for (int i = 0; i < n; i++) { rnk[ord[i]] = i; }
if (rnk[u] < rnk[v]) { dag[u].push_back(v); }
// c3
for (int x = 0; x < n; x++) {
  for (auto y : dag[x]) { vis[y] = 1; }
  for (auto y : dag[x]) { for (auto z : dag[y]) { ans += vis[z]; }}
  for (auto y : dag[x]) { vis[y] = 0; }
}
// c4
for (int x = 0; x < n; x++) {
  for (auto y : dag[x]) { for (auto z : adj[y]) { if (rnk[z] > rnk[x]) { ans += vis[z]++; }}}
  for (auto y : dag[x]) { for (auto z : adj[y]) { if (rnk[z] > rnk[x]) { vis[z]--; }}}
}
