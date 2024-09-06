struct Treap {
    Treap *lc = nullptr, *rc = nullptr;
    int sz = 1;
    unsigned w = rng();
    i64 m = 0, b = 0, val = 0;
};
int size(Treap *t) {
    return t == nullptr ? 0 : t->sz;
}
void apply(Treap *t, i64 m, i64 b) {
    t->m += m;
    t->b += b;
    t->val += m * size(t->lc) + b;
}
void pull(Treap *t) {
    t->sz = size(t->lc) + size(t->rc) + 1;
}
void push(Treap *t) {
    if (t->lc != nullptr) {
        apply(t->lc, t->m, t->b);
    }
    if (t->rc != nullptr) {
        apply(t->rc, t->m, t->b + t->m * (size(t->lc) + 1));
    }
    t->m = t->b = 0;
}
pair<Treap*, Treap*> split(Treap *t, int s) {
    if (t == nullptr) { return {t, t}; }
    push(t);
    Treap *a, *b;
    if (s <= size(t->lc)) {
        b = t;
        tie(a, b->lc) = split(t->lc, s);
    } else {
        a = t;
        tie(a->rc, b) = split(t->rc, s - size(t->lc) - 1);
    }
    pull(t);
    return {a, b};
}
Treap* merge(Treap *t1, Treap *t2) {
    if (t1 == nullptr) { return t2; }
    if (t2 == nullptr) { return t1; }
    push(t1), push(t2);
    if (t1->w > t2->w) {
        t1->rc = merge(t1->rc, t2);
        pull(t1);
        return t1;
    } else {
        t2->lc = merge(t1, t2->lc);
        pull(t2);
        return t2;
    }
}
int rnk(Treap *t, i64 val) {
    int res = 0;
    while (t != nullptr) {
        push(t);
        if (val <= t->val) {
            res += size(t->lc) + 1;
            t = t->rc;
        } else {
            t = t->lc;
        }
    }
    return res;
}
Treap* join(Treap *t1, Treap *t2) {
    if (size(t1) > size(t2)) {
        swap(t1, t2);
    }
    Treap *t = nullptr;
    while (t1 != nullptr) {
        auto [u1, v1] = split(t1, 1);
        t1 = v1;
        int r = rnk(t2, u1->val);
        if (r > 0) {
            auto [u2, v2] = split(t2, r);
            t = merge(t, u2);
            t2 = v2;
        }
        t = merge(t, u1);
    }
    t = merge(t, t2);
    return t;
}
