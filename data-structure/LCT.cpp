struct Splay {
  array<Splay*, 2> ch = {nullptr, nullptr};
  Splay* fa = nullptr;
  int sz = 1;
  bool rev = false;
  Splay() {}
  void applyRev(bool x) {
    if (x) {
      swap(ch[0], ch[1]);
      rev ^= 1;
    }
  }
  void push() {
    for (auto k : ch) {
      if (k) {
        k->applyRev(rev);
      }
    }
    rev = false;
  }
  void pull() {
    sz = 1;
    for (auto k : ch) {
      if (k) {
      }
    } 
  }
  int relation() { return this == fa->ch[1]; }
  bool isRoot() { return !fa || fa->ch[0] != this && fa->ch[1] != this; }
  void rotate() {
    Splay *p = fa;
    bool x = !relation();
    p->ch[!x] = ch[x];
    if (ch[x]) { ch[x]->fa = p; }
    fa = p->fa;
    if (!p->isRoot()) { p->fa->ch[p->relation()] = this; }
    ch[x] = p;
    p->fa = this;
    p->pull();
  }
  void splay() {
    vector<Splay*> s;
    for (Splay *p = this; !p->isRoot(); p = p->fa) { s.push_back(p->fa); }
    while (!s.empty()) {
      s.back()->push();
      s.pop_back();
    }
    push();
    while (!isRoot()) {
      if (!fa->isRoot()) {
        if (relation() == fa->relation()) {
          fa->rotate();
        } else {
          rotate();
        }
      }
      rotate();
    }
    pull();
  }
  void access() {
    for (Splay *p = this, *q = nullptr; p; q = p, p = p->fa) {
      p->splay();
      p->ch[1] = q;
      p->pull();
    }
    splay();
  }
  void makeRoot() {
    access();
    applyRev(true);
  }
  Splay* findRoot() {
    access();
    Splay *p = this;
    while (p->ch[0]) { p = p->ch[0]; }
    p->splay();
    return p;
  }
  friend void split(Splay *x, Splay *y) {
    x->makeRoot();
    y->access();
  }
  // link if not connected
  friend void link(Splay *x, Splay *y) {
    x->makeRoot();
    if (y->findRoot() != x) {
      x->fa = y;
    }
  }
  // delete edge if doesn't exist
  friend void cut(Splay *x, Splay *y) {
    split(x, y);
    if (x->fa == y && !x->ch[1]) {
      x->fa = y->ch[0] = nullptr;
      x->pull();
    }
  }
  bool connected(Splay *x, Splay *y) {
    return x->findRoot() == y->findRoot();
  }
};
