struct Treap {
  array<Treap*, 2> ch = {nullptr, nullptr};
  Treap *fa = nullptr;
  int x, P;
  int sz = 1;
  bool rev = false;

  Treap(int x = 0) : x(x), P(rng()) {}
  friend int size(Treap* t) {
    return t ? t->sz : 0;
  }
  void apply() {
    rev ^= 1;
  }
  void push() {
    if (rev) {
      swap(ch[0], ch[1]);
      for (auto k : ch) {
        if (k) {
          k->apply();
        }
      }
      rev = false;
    }
  }
  void pull() {
    sz = 1;
    for (auto *k : ch) {
      if (k) {
        sz += k->sz;
        k->fa = this;
      }
    }
  }
};
 
Treap* merge(Treap *l, Treap *r) {
  if (!l) { return r; }
  if (!r) { return l; }
  if (l->P > r->P) {
    l->push();
    l->ch[1] = merge(l->ch[1], r);
    l->pull();
    return l;
  } else {
    r->push();
    r->ch[0] = merge(l, r->ch[0]);
    r->pull();
    return r;
  }
}

pair<Treap*, Treap*> splitSize(Treap *t, int left) {
  if (t) { t->fa = nullptr; }
  if (size(t) <= left) { return {t, nullptr}; } 
  t->push();
  Treap* a;
  Treap* b;
  int sl = size(t->ch[0]) + 1;
  if (sl <= left) {
    a = t;
    tie(a->ch[1], b) = splitSize(t->ch[1], left - sl);
  } else {
    b = t;
    tie(a, b->ch[0]) = splitSize(t->ch[0], left);
  }
  t->pull();
  return {a, b};
}
