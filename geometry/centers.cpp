// radius: (a + b + c) * r / 2 = A
P<dbl> inCenter(P<dbl> a, P<dbl> b, P<dbl> c) {
    dbl la = length(b - c), lb = length(c - a), lc = length(a - b);
    return (a * la + b * lb + c * lc) / (la + lb + lc);
}

P<dbl> circumCenter(P<dbl> a, P<dbl> b, P<dbl> c) {
    P<dbl> ba = b - a, ca = c - a;
    dbl db = square(ba), dc = square(ca), d = 2 * cross(ba, ca);
    return a - P<dbl>(ba.y * dc - ca.y * db, ca.x * db - ba.x * dc) / d;
}

// radius: point to line
P<dbl> orthoCenter(P<dbl> a, P<dbl> b, P<dbl> c) {
    P ba = b - a, ca = c - a, bc = b - c;
    dbl Y = ba.y * ca.y * bc.y;
    dbl A = ca.x * ba.y - ba.x * ca.y;
    dbl x0 = (Y + ca.x * ba.y * b.x - ba.x * ca.y * c.x) / A;
    dbl y0 = -ba.x * (x0 - c.x) / ba.y + ca.y;
    return {x0, y0};
}