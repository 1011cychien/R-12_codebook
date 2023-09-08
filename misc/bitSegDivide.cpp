auto get = [&](i64 l, i64 r) {
    vector<pair<i64, i64>> res;
    if (l == 0) {
        i64 high = 1;
        while (i128(high) * 2 <= r) {
            high *= 2;
        }
        res.emplace_back(0, high - 1);
        l = high;
    }
    while (l <= r) {
        i64 nxt = l + lowbit(l) - 1;
        if (nxt > r) {
            for (int b = __builtin_ctzll(l) - 1; b >= 0; --b){
                if (l + (1ll << b) - 1 <= r){
                    res.emplace_back(l, l + (1ll << b) - 1);
                    l += 1ll << b;
                }
            }
            break;
        }
        else {
            res.emplace_back(l, nxt);
            l = nxt + 1;
        }
    }
    return res;
};
vector<pair<i64, i64>> all;
for (auto [l1, r1] : sega) {
    for (auto [l2, r2] : segb) {
        i64 length_1 = __lg(r1 - l1 + 1), length_2 = __lg(r2 - l2 + 1);

        i64 length = max(length_1, length_2);

        i64 common_prefix = ((l1 ^ l2) >> length) << length;

        i64 L = common_prefix, R = common_prefix + (1ll << length) - 1;
        all.emplace_back(L, R);
    }
}