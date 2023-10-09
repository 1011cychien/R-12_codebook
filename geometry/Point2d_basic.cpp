using Real = double; // modify these if needed
constexpr Real eps = 1e-9;
int sign(T x) { return (x > 0) - (x < 0); }
int sign(Real x) { return (x > eps) - (x < -eps); }
int cmp(T a, T b) { return sign(a - b); }
struct P {
    T x = 0, y = 0;
    P(T x = 0, T y = 0) : x(x), y(y) {}
    -, +*/, ==!=<, - (unary)
};
struct L {
    P<T> a, b;
    L(P<T> a = {}, P<T> b = {}) : a(a), b(b) {}
};
T dot(P<T> a, P<T> b) { return a.x * b.x + a.y * b.y; }
T square(P<T> a) { return dot(a, a); }
Real length(P<T> a) { return sqrtl(square(a)); }
Real dist(P<T> a, P<T> b) { return length(a - b); }
T cross(P<T> a, P<T> b) { return a.x * b.y - a.y * b.x; }
T cross(P<T> p, P<T> a, P<T> b) { return cross(a - p, b - p); }
P<Real> normal(P<T> a) {
    Real len = length(a);
    return P<Real>(a.x / len, a.y / len);
}
bool up(P<T> a) { return sign(a.y) > 0 || sign(a.y) == 0 && sign(a.x) > 0; }
// 3 colinear? please remember to remove (0, 0)
bool polar(P<T> a, P<T> b) {
    bool ua = up(a), ub = up(b);
    return ua != ub ? ua : sign(cross(a, b)) == 1;
}
bool sameDirection(P<T> a, P<T> b) { return sign(cross(a, b)) == 0 && sign(dot(a, b)) == 1; }
// 1/0/1 if on a->b's left/ /right
int side(P<T> p, P<T> a, P<T> b) { return sign(cross(p, a, b)); }
int side(P<T> p, L<T> l) { return side(p, l.a, l.b); }
P<T> rotate90(P<T> p) { return {-p.y, p.x}; }
P<Real> rotate(P<Real> p, Real ang) { return {p.x * cos(ang) - p.y * sin(ang), p.x * sin(ang) + p.y * cos(ang)}; }
Real angle(P<T> p) { return atan2(p.y, p.x); }
P<T> direction(L<T> l) { return l.b - l.a; }
bool sameDirection(L<T> l1, L<T> l2) { return sameDirection(direction(l1), direction(l2)); }
P<Real> projection(P<Real> p, L<Real> l) {
    auto d = direction(l);
    return l.a + d * (dot(p - l.a, d) / square(d));
}
P<Real> reflection(P<Real> p, L<Real> l) { return projection(p, l) * 2 - p; }
Real pointToLineDist(P<Real> p, L<Real> l) { return dist(p, projection(p, l)); }
// better use integers if you don't need exact coordinate
// l <= r is not explicitly required
P<Real> lineIntersection(L<T> l1, L<T> l2) { return l1.a - direction(l1) * (Real(cross(direction(l2), l1.a - l2.a)) / cross(direction(l2), direction(l1))); }
bool between(T m, T l, T r) { return cmp(l, m) == 0 || cmp(m, r) == 0 || l < m != r < m; }
bool pointOnSeg(P<T> p, L<T> l) { return side(p, l) == 0 && between(p.x, l.a.x, l.b.x) && between(p.y, l.a.y, l.b.y); }
bool pointStrictlyOnSeg(P<T> p, L<T> l) { return side(p, l) == 0 && sign(dot(p - l.a, direction(l))) * sign(dot(p - l.b, direction(l))) < 0; }
bool overlap(T l1, T r1, T l2, T r2) {
    if (l1 > r1) { swap(l1, r1); }
    if (l2 > r2) { swap(l2, r2); }
    return cmp(r1, l2) != -1 && cmp(r2, l1) != -1;
}
bool segIntersect(L<T> l1, L<T> l2) {
    auto [p1, p2] = l1;
    auto [q1, q2] = l2;
    return overlap(p1.x, p2.x, q1.x, q2.x) && overlap(p1.y, p2.y, q1.y, q2.y) &&
            side(p1, l2) * side(p2, l2) <= 0 && 
            side(q1, l1) * side(q2, l1) <= 0;
}
// parallel intersecting is false
bool segStrictlyIntersect(L<T> l1, L<T> l2) {
    auto [p1, p2] = l1;
    auto [q1, q2] = l2;
    return side(p1, l2) * side(p2, l2) < 0 && 
           side(q1, l1) * side(q2, l1) < 0;
}
// parallel or intersect at source doesn't count
bool rayIntersect(L<T> l1, L<T> l2) {
    int x = sign(cross(l1.b - l1.a, l2.b - l2.a));
    return x == 0 ? false : side(l1.a, l2) == x && side(l2.a, l1) == -x;
}
Real pointToSegDist(P<T> p, L<T> l) {
    P<Real> q = projection(p, l);
    if (pointOnSeg(q, l)) {
        return dist(p, q);
    } else {
        return min(dist(p, l.a), dist(p, l.b));
    }
}
Real segDist(L<T> l1, L<T> l2) {
    if (segIntersect(l1, l2)) { return 0; }
	return min({pointToSegDist(l1.a, l2), pointToSegDist(l1.b, l2), 
            pointToSegDist(l2.a, l1), pointToSegDist(l2.b, l1)});
}
// 2 times area
T area(vector<P<T>> a) {
    T res = 0;
    int n = a.size();
    for (int i = 0; i < n; i++) { res += cross(a[i], a[(i + 1) % n]); }
    return res;
}
bool pointInPoly(P<T> p, vector<P<T>> a) {
    int n = a.size(), res = 0;
    for (int i = 0; i < n; i++) {
        P<T> u = a[i], v = a[(i + 1) % n];
        if (pointOnSeg(p, {u, v})) { return 1; }
        if (cmp(u.y, v.y) <= 0) { swap(u, v); }
        if (cmp(p.y, u.y) > 0 || cmp(p.y, v.y) <= 0) { continue; }
        res ^= cross(p, u, v) > 0;
    }
    return res;
}