vector<int> vis(n);
// c3
for (int x = 0; x < n; x++) {
    for (auto y : dag[x]) {
        vis[y] = 1;
    }
    for (auto y : dag[x]) {
        for (auto z : dag[y]) {
            ans += vis[z];
        }
    }
    for (auto y : dag[x]) {
        vis[y] = 0;
    }
}
// c4
for (int x = 0; x < n; x++) {
    for (auto y : dag[x]) {
        for (auto z : adj[y]) {
            if (z != x) {
                ans += vis[z]++;
            }
        }
    }
    for (auto y : dag[x]) {
        for (auto z : adj[y]) {
            if (z != x) {
                vis[z]--;
            }
        }
    }
}