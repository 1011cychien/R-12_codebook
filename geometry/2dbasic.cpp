using dbl = double; // or long double

constexpr dbl eps = 1e-9;
int sign(dbl x) {
    return (x > eps) - (x < -eps);
}
int dcmp(dbl a, dbl b) {
    return sign(a - b);
}

template <typename T>
struct P {
    T x = 0, y = 0;
    P(T x = 0, T y = 0) : x(x), y(y) {}
    friend istream &operator>>(istream &is, P &p) { return is >> p.x >> p.y; }
    friend ostream &operator<<(ostream &os, P p) { return os << '(' << p.x << ", " << p.y << ')'; }
    friend bool operator==(P a, P b) { return dcmp(a.x, b.x) == 0 && dcmp(a.y, b.y) == 0; }
    friend bool operator!=(P a, P b) { return !(a == b); }
    P operator-() { return P(-x, -y); }
    P& operator+=(P a) {
        x += a.x, y += a.y;
        return *this;
    }
    P& operator-=(P a) {
        x -= a.x, y -= a.y;
        return *this;
    }
    P& operator*=(T d) {
        x *= d, y *= d;
        return *this;
    }
    P& operator/=(T d) {
        x /= d, y /= d;
        return *this;
    }
    friend P operator+(P a, P b) { return P(a) += b; }
    friend P operator-(P a, P b) { return P(a) -= b; }
    friend P operator*(P a, T d) { return P(a) *= d; }
    friend P operator/(P a, T d) { return P(a) /= d; }
    friend bool operator<(P a, P b) {
        int sx = dcmp(a.x, b.x);
        return sx != 0 ? sx == -1 : dcmp(a.y, b.y) == -1;
    }
};

template <typename T>
struct L {
    P<T> a, b;
    L(P<T> a = {}, P<T> b = {}) : a(a), b(b) {}
};

template <typename T>
T dot(P<T> a, P<T> b) { return a.x * b.x + a.y * b.y; }

template <typename T>
T square(P<T> a) { return dot(a, a); }

template <typename T>
dbl length(P<T> a) { return sqrtl(square(a)); }

template <typename T>
T cross(P<T> a, P<T> b) { return a.x * b.y - a.y * b.x; }

template <typename T>
T cross(P<T> p, P<T> a, P<T> b) { return cross(a - p, b - p); }

template <typename T>
bool up(P<T> a) { return sign(a.y) > 0 || sign(a.y) == 0 && sign(a.x) > 0; }

// 3 colinear? please remember to remove (0, 0)
template <typename T>
bool polar(P<T> a, P<T> b) {
    bool ua = up(a), ub = up(b);
    return ua != ub ? ua : sign(cross(a, b)) == 1;
};

// 1 if on a->b's left
template <typename T>
int side(P<T> p, P<T> a, P<T> b) { return sign(cross(p, a, b)); }

template <typename T>
int side(P<T> p, L<T> l) { return side(p, l.a, l.b); }

