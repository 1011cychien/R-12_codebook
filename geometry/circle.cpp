const Real PI = acos(-1);
struct Circle {
  P<Real> o;
  Real r;
  Circle(P<Real> o = {}, Real r = 0) : o(o), r(r) {}
};
// actually counts number of tangent lines
int typeOfCircles(Circle c1, Circle c2) {
  auto [o1, r1] = c1;
  auto [o2, r2] = c2;
  Real d = dist(o1, o2);
  if (cmp(d, r1 + r2) == 1) { return 4; }
  if (cmp(d, r1 + r2) == 0) { return 3; }
  if (cmp(d, abs(r1 - r2)) == 1) { return 2; }
  if (cmp(d, abs(r1 - r2)) == 0) { return 1; }
  return 0;
}
// aligned l.a -> l.b;
vector<P<Real>> circleLineIntersection(Circle c, L<Real> l) {
  P<Real> p = projection(c.o, l);
  Real h = c.r * c.r - square(p - c.o);
  if (sign(h) < 0) { return {}; }
  P<Real> q = normal(direction(l)) * sqrtl(c.r * c.r - square(p - c.o));
  return {p - q, p + q};
}
// circles shouldn't be identical
// duplicated if only one intersection, aligned c1 counterclockwise
vector<P<Real>> circleIntersection(Circle c1, Circle c2) {
  int type = typeOfCircles(c1, c2);
  if (type == 0 || type == 4) { return {}; }
  auto [o1, r1] = c1;
  auto [o2, r2] = c2;
  Real d = clamp(dist(o1, o2), abs(r1 - r2), r1 + r2);
  Real y = (r1 * r1 + d * d - r2 * r2) / (2 * d), x = sqrtl(r1 * r1 - y * y);
  P<Real> dir = normal(o2 - o1), q1 = o1 + dir * y, q2 = rotate90(dir) * x;
  return {q1 - q2, q1 + q2};
}
// counterclockwise, on circle -> no tangent
vector<P<Real>> pointCircleTangent(P<Real> p, Circle c) {
  Real x = square(p - c.o), d = x - c.r * c.r;
  if (sign(d) <= 0) { return {}; }
  P<Real> q1 = c.o + (p - c.o) * (c.r * c.r / x), q2 = rotate90(p - c.o) * (c.r * sqrt(d) / x);
  return {q1 - q2, q1 + q2};
}
// one-point tangent lines are not returned
vector<L<Real>> externalTangent(Circle c1, Circle c2) {
  auto [o1, r1] = c1;
  auto [o2, r2] = c2;
  vector<L<Real>> res;
  if (cmp(r1, r2) == 0) {
    P dr = rotate90(normal(o2 - o1)) * r1;
    res.emplace_back(o1 + dr, o2 + dr);
    res.emplace_back(o1 - dr, o2 - dr);
  } else {
    P p = (o2 * r1 - o1 * r2) / (r1 - r2);
    auto ps = pointCircleTangent(p, c1), qs = pointCircleTangent(p, c2);
    for (int i = 0; i < int(min(ps.size(), qs.size())); i++) { res.emplace_back(ps[i], qs[i]); }
  }
  return res;
}
vector<L<Real>> internalTangent(Circle c1, Circle c2) {
  auto [o1, r1] = c1;
  auto [o2, r2] = c2;
  vector<L<Real>> res;
  P<Real> p = (o1 * r2 + o2 * r1) / (r1 + r2);
  auto ps = pointCircleTangent(p, c1), qs = pointCircleTangent(p, c2);
  for (int i = 0; i < int(min(ps.size(), qs.size())); i++) { res.emplace_back(ps[i], qs[i]); }
  return res;
}
// OAB and circle directed area
Real triangleCircleIntersectionArea(P<Real> p1, P<Real> p2, Real r) {
  auto angle = [&](P<Real> p1, P<Real> p2) { return atan2l(cross(p1, p2), dot(p1, p2)); };
  vector<P<Real>> v = circleLineIntersection(Circle(P<Real>(), r), L<Real>(p1, p2));
  if (v.empty()) { return r * r * angle(p1, p2) / 2; }
  bool b1 = cmp(square(p1), r * r) == 1, b2 = cmp(square(p2), r * r) == 1;
  if (b1 && b2) {
    if (sign(dot(p1 - v[0], p2 - v[0])) <= 0 && sign(dot(p1 - v[0], p2 - v[0])) <= 0) {
      return r * r * (angle(p1, v[0]) + angle(v[1], p2)) / 2 + cross(v[0], v[1]) / 2;
    } else {
      return r * r * angle(p1, p2) / 2;
    }
  } else if (b1) {
    return (r * r * angle(p1, v[0]) + cross(v[0], p2)) / 2;
  } else if (b2) {
    return (cross(p1, v[1]) + r * r * angle(v[1], p2)) / 2;
  } else {
    return cross(p1, p2) / 2;
  }
}
Real polyCircleIntersectionArea(const vector<P<Real>> &a, Circle c) {
  int n = a.size();
  Real ans = 0;
  for (int i = 0; i < n; i++) {
    ans += triangleCircleIntersectionArea(a[i], a[(i + 1) % n], c.r);
  }
  return ans;
}
Real circleIntersectionArea(Circle a, Circle b) {
  int t = typeOfCircles(a, b);
  if (t >= 3) {
    return 0;
  } else if (t <= 1) {
    Real r = min(a.r, b.r);
    return r * r * PI;
  }
  Real res = 0, d = dist(a.o, b.o);
  for (int i = 0; i < 2; ++i) {
    Real alpha = acos((b.r * b.r + d * d - a.r * a.r) / (2 * b.r * d));
    Real s = alpha * b.r * b.r;
    Real t = b.r * b.r * sin(alpha) * cos(alpha);
    res += s - t;
    swap(a, b);
  }
  return res;
}
