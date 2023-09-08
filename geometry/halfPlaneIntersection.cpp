vector<P<Real>> halfPlaneIntersection(vector<L<Real>> a) {
    sort(a.begin(), a.end(), [&](auto l1, auto l2) { 
        if (sameDirection(l1, l2)) {
            return side(l1.a, l2) > 0;
        } else {
            return polar(direction(l1), direction(l2));
        }
    });
    deque<L<Real>> dq;
    auto check = [&](L<Real> l, L<Real> l1, L<Real> l2) { return side(lineIntersection(l1, l2), l) > 0; };
    for (int i = 0; i < int(a.size()); i++) {
        if (i > 0 && sameDirection(a[i], a[i - 1])) { continue; }
        while (int(dq.size()) > 1 && !check(a[i], dq.end()[-2], dq.back())) { dq.pop_back(); }
        while (int(dq.size()) > 1 && !check(a[i], dq[1], dq[0])) { dq.pop_front(); }
        dq.push_back(a[i]);
    }
    while (int(dq.size()) > 2 && !check(dq[0], dq.end()[-2], dq.back())) { dq.pop_back(); }
    while (int(dq.size()) > 2 && !check(dq.back(), dq[1], dq[0])) { dq.pop_front(); }
    vector<P<Real>> res;
    dq.push_back(dq[0]);
    for (int i = 0; i + 1 < int(dq.size()); i++) { res.push_back(lineIntersection(dq[i], dq[i + 1])); }
    return res;
}