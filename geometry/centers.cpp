// radius: (a + b + c) * r / 2 = A or pointToLineDist
P<Real> inCenter(P<Real> a, P<Real> b, P<Real> c) {
    Real la = length(b - c), lb = length(c - a), lc = length(a - b);
    return (a * la + b * lb + c * lc) / (la + lb + lc);
}
// used in min enclosing circle
P<Real> circumCenter(P<Real> a, P<Real> b, P<Real> c) {
    P<Real> ba = b - a, ca = c - a;
    Real db = square(ba), dc = square(ca), d = 2 * cross(ba, ca);
    return a - P<Real>(ba.y * dc - ca.y * db, ca.x * db - ba.x * dc) / d;
}
P<Real> orthoCenter(P<Real> a, P<Real> b, P<Real> c) {
    L<Real> u(c, P<Real>(c.x - a.y + b.y, c.y + a.x - b.x));
    L<Real> v(b, P<Real>(b.x - a.y + c.y, b.y + a.x - c.x));
    return lineIntersection(u, v);
}
