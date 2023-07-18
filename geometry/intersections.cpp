// better use integers when dealing with intersections
// l <= r is not explicitly required

template <typename T>
bool between(T l, T m, T r) {
    return dcmp(l, m) == 0 || dcmp(m, r) == 0 || (l < m != r < m);
}

// template <typename T>
// bool pointOnSegment() {

// }

// template <typename T>
// bool pointStrictlyOnSegment() {

// }

template <typename T>
bool overlap(T l1, T r1, T l2, T r2) {
    if (l1 > r1) { swap(l1, r1); }
    if (l2 > r2) { swap(l2, r2); }
    return dcmp(r1, l2) != -1 && dcmp(r2, l1) != -1;
}

template <typename T>
bool segmentIntersect(L<T> l1, L<T> l2) {
    auto [p1, p2] = l1;
    auto [q1, q2] = l2;
    return overlap(p1.x, p2.x, q1.x, q2.x) && overlap(p1.y, p2.y, q1.y, q2.y) &&
           side(p1, p2, q1) * side(p1, p2, q2) <= 0 && 
           side(q1, q2, p1) * side(q1, q2, p2) <= 0;
}

// parallel intersecting is false
template <typename T>
bool segmentStrictlyIntersect(L<T> l1, L<T> l2) {
    auto [p1, p2] = l1;
    auto [q1, q2] = l2;
    return side(p1, p2, q1) * side(p1, p2, q2) < 0 && 
           side(q1, q2, p1) * side(q1, q2, p2) < 0;
}

// parallel or intersect at source doesn't count
template <typename T>
bool rayIntersect(L<T> l1, L<T> l2) {
  int x = sign(cross(l1.b - l1.a, l2.b - l2.a));
  return x == 0 ? false : side(l1.a, l2) == x && side(l2.a, l1) == -x;
}