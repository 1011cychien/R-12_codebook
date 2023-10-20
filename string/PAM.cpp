// cnt : occurrences, (dfs fail tree)
// num : number of pal ending here
struct PAM {
  static constexpr int A = 26;
  struct Node {
    int len = 0, link = 0, cnt = 0, num = 0;
    array<int, A> nxt{};
    Node() {}
  };
  vector<Node> t;
  int suf = 1;
  string s;
  PAM() : t(2) { t[0].len = -1; }
  int size() { return t.size(); }
  Node& operator[](int i) { return t[i]; }
  int newNode() {
    t.emplace_back();
    return t.size() - 1;
  }
  bool add(int c, char offset = 'a') {
    int pos = s.size();
    s += c + offset;
    int cur = suf, curlen = 0;
    while (true) {
      curlen = t[cur].len;
      if (pos - 1 - curlen >= 0 && s[pos - 1 - curlen] == s[pos]) { break; }
      cur = t[cur].link;
    }
    if (t[cur].nxt[c]) {
      suf = t[cur].nxt[c];
      t[suf].cnt++;
      return false;
    }
    suf = newNode();
    t[suf].len = t[cur].len + 2;
    t[suf].cnt = t[suf].num = 1;
    t[cur].nxt[c] = suf;
    if (t[suf].len == 1) {
      t[suf].link = 1;
      return true;
    }
    while (true) {
      cur = t[cur].link;
      curlen = t[cur].len;
      if (pos - 1 - curlen >= 0 && s[pos - 1 - curlen] == s[pos]) {
        t[suf].link = t[cur].nxt[c];
        break;
      }
    }
    t[suf].num += t[t[suf].link].num;
    return true;
  }
};
