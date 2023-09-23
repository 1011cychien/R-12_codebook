// find \sum a_(i- j)c_j = 0 for d <= i
template <typename T>
vector<T> berlekampMassey(const vector<T> &a) {
    vector<T> c(1, 1), oldC(1);
    int oldI = -1;
    T oldD = 1;
    for (int i = 0; i < int(a.size()); i++) {
        T d = 0;
        for (int j = 0; j < int(c.size()); j++) { d += c[j] * a[i - j]; }
        if (d == 0) { continue; }
        T mul = d / oldD;
        vector<T> nc = c;
        nc.resize(max(int(c.size()), i - oldI + int(oldC.size())));
        for (int j = 0; j < int(oldC.size()); j++) { nc[j + i - oldI] -= oldC[j] * mul; }
        if (i - int(c.size()) > oldI - int(oldC.size())) {
            oldI = i; 
            oldD = d;
            swap(oldC, c);
        }
        swap(c, nc);
    }
    return c;
}