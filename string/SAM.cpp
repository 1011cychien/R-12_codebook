struct SAM {
  static constexpr int A = 26;
  struct Node {
    int len = 0, link = -1, cnt = 0;
    array<int, A> nxt;
    Node() { nxt.fill(-1); }
  };
  vector<Node> t;
  SAM() : t(1) {}
  int size() { return t.size(); }
  Node& operator[](int i) { return t[i]; }
  int newNode() {
    t.emplace_back();
    return t.size() - 1;
  }
  int extend(int p, int c) {
    int cur = newNode();
    t[cur].len = t[p].len + 1;
    t[cur].cnt = 1;
    while (p != -1 && t[p].nxt[c] == -1) {
      t[p].nxt[c] = cur;
      p = t[p].link;
    }
    if (p == -1) {
      t[cur].link = 0;
    } else {
      int q = t[p].nxt[c];
      if (t[p].len + 1 == t[q].len) {
        t[cur].link = q;
      } else {
        int clone = newNode();
        t[clone].len = t[p].len + 1;
        t[clone].link = t[q].link;
        t[clone].nxt = t[q].nxt;
        while (p != -1 && t[p].nxt[c] == q) {
          t[p].nxt[c] = clone;
          p = t[p].link;
        }
        t[q].link = t[cur].link = clone;
      }
    }
    return cur;
  }
};
