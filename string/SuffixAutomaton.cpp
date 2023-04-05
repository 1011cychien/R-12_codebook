constexpr int K = 26;
struct Node{
    int len = 0, link = -1, cnt = 0;
    array<int, K> nxt;
    Node() { nxt.fill(-1); }
};
vector<Node> sam(1);
auto extend = [&](int c) {
    static int last = 0;
    int p = last, cur = sam.size();
    sam.emplace_back();
    sam[cur].len = sam[p].len + 1;
    sam[cur].cnt = 1;
    while (p != -1 && sam[p].nxt[c] == -1) {
        sam[p].nxt[c] = cur;
        p = sam[p].link;
    }
    if (p == -1) {
        sam[cur].link = 0;
    } else {
        int q = sam[p].nxt[c];
        if (sam[p].len + 1 == sam[q].len) {
            sam[cur].link = q;
        } else {
            int clone = sam.size();
            sam.emplace_back();
            sam[clone].len = sam[p].len + 1;
            sam[clone].link = sam[q].link;
            sam[clone].nxt = sam[q].nxt;
            while (p != -1 && sam[p].nxt[c] == q) {
                sam[p].nxt[c] = clone;
                p = sam[p].link;
            }
            sam[q].link = sam[cur].link = clone;
        }
    }
    last = cur;
};
for (auto ch : s) {
    extend(ch - 'a');
}
int N = sam.size();
vector<vector<int>> g(N);
for (int i = 1; i < N; i++) {
    g[sam[i].link].push_back(i);
}