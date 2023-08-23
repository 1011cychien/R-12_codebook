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
    P ba = b - a, ca = c - a, bc = b - c;
    Real Y = ba.y * ca.y * bc.y;
    Real A = ca.x * ba.y - ba.x * ca.y;
    Real x0 = (Y + ca.x * ba.y * b.x - ba.x * ca.y * c.x) / A;
    Real y0 = -ba.x * (x0 - c.x) / ba.y + ca.y;
    return {x0, y0};
}