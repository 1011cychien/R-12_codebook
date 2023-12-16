bool disjunct(Circle &a, Circle &b, bool strict) { return sign(dist(a.o, b.o) - a.r - b.r) > (strict ? 0 : -1); }
bool contain(Circle &a, Circle &b, bool strict) { return sign(a.r - b.r - dist(a.o, b.o)) > (strict ? 0 : -1); }
const Real PI = acos(-1);

Real circleUnion(vector<Circle> &c) {
    int n = c.size();
    Real area = 0;
    vector overlap(n, vector<int>(n)), g(overlap);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            overlap[i][j] = i != j && (cmp(c[i].r, c[j].r) > 0 || cmp(c[i].r, c[j].r) == 0 && i < j) && contain(c[i], c[j], false);
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            g[i][j] = i != j && !(overlap[i][j] || overlap[j][i] || disjunct(c[i], c[j], false));
        }
    }
    for (int i = 0; i < n; i++) {
        vector<Teve> events;
        int cnt = 1;
        for (int j = 0; j < n; j++) {
            cnt += overlap[j][i];
        }
        for (int j = 0; j < n; j++) {
            if (g[i][j]) {
                auto pts = circleIntersection(c[i], c[j]);
                auto a = pts[1], b = pts[0];
                Real A = angle(a - c[i].o), B = angle(b - c[i].o);
                events.emplace_back(b, B, 1);
                events.emplace_back(a, A, -1);
                if (B > A) {
                    cnt++;
                }
            }
        }
        if (events.empty()) {
            if (cnt == 1) {
                area += PI * c[i].r * c[i].r;
            }
        } else {
            sort(events.begin(), events.end());
            events.emplace_back(events[0]);
            events.back().ang += 2 * PI;
            for (int j = 0; j + 1 < int(events.size()); j++) {
                cnt += events[j].add;
                if (cnt == 1) { area += cross(events[j].p, events[j + 1].p) * 0.5; }
                Real t = events[j + 1].ang - events[j].ang;
                if (cnt == 1) { area += (t - sin(t)) * c[i].r * c[i].r * 0.5; }
            }
        }
    }
    return area;
}