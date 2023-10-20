const P<i64> pINF = P<i64>(1e18, 1e18);
using i128 = __int128_t;
struct Quad {
  P<i64> origin;
  Quad *rot = nullptr, *onext = nullptr;
  bool used = false;
  Quad* rev() const { return rot->rot; }
  Quad* lnext() const { return rot->rev()->onext->rot; }
  Quad* oprev() const { return rot->onext->rot; }
  P<i64> dest() const { return rev()->origin; }
};
Quad* makeEdge(P<i64> from, P<i64> to) {
  Quad *e1 = new Quad, *e2 = new Quad, *e3 = new Quad, *e4 = new Quad;
  e1->origin = from;
  e2->origin = to;
  e3->origin = e4->origin = pINF;
  e1->rot = e3;
  e2->rot = e4;
  e3->rot = e2;
  e4->rot = e1;
  e1->onext = e1;
  e2->onext = e2;
  e3->onext = e4;
  e4->onext = e3;
  return e1;
}
void splice(Quad *a, Quad *b) {
  swap(a->onext->rot->onext, b->onext->rot->onext);
  swap(a->onext, b->onext);
}
void delEdge(Quad *e) {
  splice(e, e->oprev());
  splice(e->rev(), e->rev()->oprev());
  delete e->rev()->rot;
  delete e->rev();
  delete e->rot;
  delete e;
}
Quad *connect(Quad *a, Quad *b) {
  Quad *e = makeEdge(a->dest(), b->origin);
  splice(e, a->lnext());
  splice(e->rev(), b);
  return e;
}
bool onLeft(P<i64> p, Quad *e) { return side(p, e->origin, e->dest()) > 0; }
bool onRight(P<i64> p, Quad *e) { return side(p, e->origin, e->dest()) < 0; }
template <class T>
T det3(T a1, T a2, T a3, T b1, T b2, T b3, T c1, T c2, T c3) {
  return a1 * (b2 * c3 - c2 * b3) - a2 * (b1 * c3 - c1 * b3) + a3 * (b1 * c2 - c1 * b2);
}
bool inCircle(P<i64> a, P<i64> b, P<i64> c, P<i64> d) {
  auto f = [&](P<i64> a, P<i64> b, P<i64> c) {
    return det3<i128>(a.x, a.y, square(a), b.x, b.y, square(b), c.x, c.y, square(c));
  };
  i128 det = f(a, c, d) + f(a, b, c) - f(b, c, d) - f(a, b, d);
  return det > 0;
}
pair<Quad*, Quad*> build(int l, int r, vector<P<i64>> &p) {
  if (r - l == 2) {
    Quad *res = makeEdge(p[l], p[l + 1]);
    return pair(res, res->rev());
  } else if (r - l == 3) {
    Quad *a = makeEdge(p[l], p[l + 1]), *b = makeEdge(p[l + 1], p[l + 2]);
    splice(a->rev(), b);
    int sg = sign(cross(p[l], p[l + 1], p[l + 2]));
    if (sg == 0) { return pair(a, b->rev()); }
    Quad *c = connect(b, a);
    if (sg == 1) {
      return pair(a, b->rev());
    } else {
      return pair(c->rev(), c);
    }
  }
  int m = l + r >> 1;
  auto [ldo, ldi] = build(l, m, p);
  auto [rdi, rdo] = build(m, r, p);
  while (true) {
    if (onLeft(rdi->origin, ldi)) {
      ldi = ldi->lnext();
      continue;
    }
    if (onRight(ldi->origin, rdi)) {
      rdi = rdi->rev()->onext;
      continue;
    }
    break;
  }
  Quad *basel = connect(rdi->rev(), ldi);
  auto valid = [&](Quad *e) { return onRight(e->dest(), basel); };
  if (ldi->origin == ldo->origin) { ldo = basel->rev(); }
  if (rdi->origin == rdo->origin) { rdo = basel; }
  while (true) {
    Quad *lcand = basel->rev()->onext;
    if (valid(lcand)) {
      while (inCircle(basel->dest(), basel->origin, lcand->dest(), lcand->onext->dest())) {
        Quad *t = lcand->onext;
        delEdge(lcand);
        lcand = t;
      }
    }
    Quad *rcand = basel->oprev();
    if (valid(rcand)) {
      while (inCircle(basel->dest(), basel->origin, rcand->dest(), rcand->oprev()->dest())) {
        Quad *t = rcand->oprev();
        delEdge(rcand);
        rcand = t;
      }
    }
    if (!valid(lcand) && !valid(rcand)) { break; }
    if (!valid(lcand) || valid(rcand) && inCircle(lcand->dest(), lcand->origin, rcand->origin, rcand->dest())) {
      basel = connect(rcand, basel->rev());
    } else {
      basel = connect(basel->rev(), lcand->rev());
    }
  }
  return pair(ldo, rdo);
}
vector<array<P<i64>, 3>> delaunay(vector<P<i64>> p) {
  sort(p.begin(), p.end());
  auto res = build(0, p.size(), p);
  Quad *e = res.first;
  vector<Quad*> edges = {e};
  while (sign(cross(e->onext->dest(), e->dest(), e->origin)) == -1) { e = e->onext; }
  auto add = [&]() {
    Quad *cur = e;
    do {
      cur->used = true;
      p.push_back(cur->origin);
      edges.push_back(cur->rev());
      cur = cur->lnext();
    } while (cur != e);
  };
  add();
  p.clear();
  int i = 0; 
  while (i < int(edges.size())) { if (!(e = edges[i++])->used) { add(); }}
  vector<array<P<i64>, 3>> ans(p.size() / 3);
  for (int i = 0; i < int(p.size()); i++) { ans[i / 3][i % 3] = p[i]; }
  return ans;
}
