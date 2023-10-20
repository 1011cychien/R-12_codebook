const int N = 1021;
struct CircleCover {
  int C; 
  Cir c[N];
  bool g[N][N], overlap[N][N];
  // Area[i] : area covered by at least i circles
  double Area[ N ];
  void init(int _C){ C = _C;}
  struct Teve {
    pdd p; double ang; int add;
    Teve() {}
    Teve(pdd _a, double _b, int _c):p(_a), ang(_b), add(_c){}
    bool operator<(const Teve &a)const
    {return ang < a.ang;}
  }eve[N * 2];
  // strict: x = 0, otherwise x = -1
  bool disjuct(Cir &a, Cir &b, int x)
  {return sign(abs(a.O - b.O) - a.R - b.R) > x;}
  bool contain(Cir &a, Cir &b, int x)
  {return sign(a.R - b.R - abs(a.O - b.O)) > x;}
  bool contain(int i, int j) {
    /* c[j] is non-strictly in c[i]. */
    return (sign(c[i].R - c[j].R) > 0 || (sign(c[i].R - c[j].R) == 0 && i < j)) && contain(c[i], c[j], -1);
  }
  void solve(){
    fill_n(Area, C + 2, 0);
    for(int i = 0; i < C; ++i)
      for(int j = 0; j < C; ++j)
        overlap[i][j] = contain(i, j);
    for(int i = 0; i < C; ++i)
      for(int j = 0; j < C; ++j) 
        g[i][j] = !(overlap[i][j] || overlap[j][i] ||
            disjuct(c[i], c[j], -1));
    for(int i = 0; i < C; ++i){
      int E = 0, cnt = 1;
      for(int j = 0; j < C; ++j)
        if(j != i && overlap[j][i])
          ++cnt;
      for(int j = 0; j < C; ++j)
        if(i != j && g[i][j]) {
          pdd aa, bb;
          CCinter(c[i], c[j], aa, bb);
          double A = atan2(aa.Y - c[i].O.Y, aa.X - c[i].O.X);
          double B = atan2(bb.Y - c[i].O.Y, bb.X - c[i].O.X);
          eve[E++] = Teve(bb, B, 1), eve[E++] = Teve(aa, A, -1);
          if(B > A) ++cnt;
        }
      if(E == 0) Area[cnt] += pi * c[i].R * c[i].R;
      else{
        sort(eve, eve + E);
        eve[E] = eve[0];
        for(int j = 0; j < E; ++j){
          cnt += eve[j].add; 
          Area[cnt] += cross(eve[j].p, eve[j + 1].p) * .5;
          double theta = eve[j + 1].ang - eve[j].ang;
          if (theta < 0) theta += 2. * pi;
          Area[cnt] += (theta - sin(theta)) * c[i].R * c[i].R * .5;
        }
      }
    }
  }
};
// p, q is convex
double TwoConvexHullMinDist(Point P[], Point Q[], int n, int m) {
  int YMinP = 0, YMaxQ = 0;
  double tmp, ans = 999999999;
  for (i = 0; i < n; ++i) if(P[i].y < P[YMinP].y) YMinP = i;
  for (i = 0; i < m; ++i) if(Q[i].y > Q[YMaxQ].y) YMaxQ = i;
  P[n] = P[0], Q[m] = Q[0];
  for (int i = 0; i < n; ++i) {
    while (tmp = Cross(Q[YMaxQ + 1] - P[YMinP + 1], P[YMinP] - P[YMinP + 1]) > Cross(Q[YMaxQ] - P[YMinP + 1], P[YMinP] - P[YMinP + 1])) YMaxQ = (YMaxQ + 1) % m;
    if (tmp < 0) ans = min(ans, PointToSegDist(P[YMinP], P[YMinP + 1], Q[YMaxQ]));
    else ans = min(ans, TwoSegMinDist(P[YMinP], P[YMinP + 1], Q[YMaxQ], Q[YMaxQ + 1]));
    YMinP = (YMinP + 1) % n;
  }
  return ans;
}
template <typename F, typename C> class MCMF {
  static constexpr F INF_F = numeric_limits<F>::max();
  static constexpr C INF_C = numeric_limits<C>::max();
  vector<tuple<int, int, F, C>> es;
  vector<vector<int>> g;
  vector<F> f;
  vector<C> d;
  vector<int> pre, inq;
  void spfa(int s) {
    fill(inq.begin(), inq.end(), 0);
    fill(d.begin(), d.end(), INF_C);
    fill(pre.begin(), pre.end(), -1);

    queue<int> q;
    d[s] = 0;
    q.push(s);
    while (!q.empty()) {
      int u = q.front();
      inq[u] = false;
      q.pop();
      for (int j : g[u]) {
        int to = get<1>(es[j]);
        C w = get<3>(es[j]);
        if (f[j] == 0 || d[to] <= d[u] + w)
          continue;
        d[to] = d[u] + w;
        pre[to] = j;
        if (!inq[to]) {
          inq[to] = true;
          q.push(to);
        }
      }
    }
  }
public:
  MCMF(int n) : g(n), pre(n), inq(n) {}
  void add_edge(int s, int t, F c, C w) {
    g[s].push_back(es.size());
    es.emplace_back(s, t, c, w);
    g[t].push_back(es.size());
    es.emplace_back(t, s, 0, -w);
  }
  pair<F, C> solve(int s, int t, C mx = INF_C / INF_F) {
    add_edge(t, s, INF_F, -mx);
    f.resize(es.size()), d.resize(es.size());
    for (F I = INF_F ^ (INF_F / 2); I; I >>= 1) {
      for (auto &fi : f)
        fi *= 2;
      for (size_t i = 0; i < f.size(); i += 2) {
        auto [u, v, c, w] = es[i];
        if ((c & I) == 0)
          continue;
        if (f[i]) {
          f[i] += 1;
          continue;
        }
        spfa(v);
        if (d[u] == INF_C || d[u] + w >= 0) {
          f[i] += 1;
          continue;
        }
        f[i + 1] += 1;
        while (u != v) {
          int x = pre[u];
          f[x] -= 1;
          f[x ^ 1] += 1;
          u = get<0>(es[x]);
        }
      }
    }

    C w = 0;
    for (size_t i = 1; i + 2 < f.size(); i += 2)
      w -= f[i] * get<3>(es[i]);
    return {f.back(), w};
  }
};
  auto [f, c] = mcmf.solve(s, t, 1e12);
  cout << f << ' ' << c << '\n';

void MoAlgoOnTree() {
  Dfs(0, -1);
  vector<int> euler(tk);
  for (int i = 0; i < n; ++i) {
    euler[tin[i]] = i;
    euler[tout[i]] = i;
  }
  vector<int> l(q), r(q), qr(q), sp(q, -1);
  for (int i = 0; i < q; ++i) {
    if (tin[u[i]] > tin[v[i]]) swap(u[i], v[i]);
    int z = GetLCA(u[i], v[i]);
    sp[i] = z[i];
    if (z == u) l[i] = tin[u[i]], r[i] = tin[v[i]];
    else l[i] = tout[u[i]], r[i] = tin[v[i]];
    qr[i] = i;
  }
  sort(qr.begin(), qr.end(), [&](int i, int j) {
      if (l[i] / kB == l[j] / kB) return r[i] < r[j];
      return l[i] / kB < l[j] / kB;
      });
  vector<bool> used(n);
  // Add(v): add/remove v to/from the path based on used[v]
  for (int i = 0, tl = 0, tr = -1; i < q; ++i) {
    while (tl < l[qr[i]]) Add(euler[tl++]);
    while (tl > l[qr[i]]) Add(euler[--tl]);
    while (tr > r[qr[i]]) Add(euler[tr--]);
    while (tr < r[qr[i]]) Add(euler[++tr]);
    // add/remove LCA(u, v) if necessary
  }
}

for (int l = 0, r = -1; auto [ql, qr, i] : qs) {
    if (ql / B == qr / B) {
        for (int j = ql; j <= qr; j++) {
            cntSmall[a[j]]++;
            ans[i] = max(ans[i], 1LL * b[a[j]] * cntSmall[a[j]]);
        }
        for (int j = ql; j <= qr; j++) {
            cntSmall[a[j]]--;
        }
        continue;
    }
    if (int block = ql / B; block != lst) {
        int x = min((block + 1) * B, n);
        while (r + 1 < x) { add(++r); }
        while (r >= x) { del(r--); }
        while (l < x) { del(l++); }
        mx = 0;
        lst = block;
    }
    while (r < qr) { add(++r); }
    i64 tmpMx = mx;
    int tmpL = l;
    while (l > ql) { add(--l); }
    ans[i] = mx;
    mx = tmpMx;
    while (l < tmpL) { del(l++); }
}

typedef pair<ll,int> T;
typedef struct heap* ph;
struct heap { // min heap
	ph l = NULL, r = NULL;
	int s = 0; T v; // s: path to leaf
	heap(T _v):v(_v) {}
};
ph meld(ph p, ph q) {
	if (!p || !q) return p?:q;
	if (p->v > q->v) swap(p,q);
	ph P = new heap(*p); P->r = meld(P->r,q);
	if (!P->l || P->l->s < P->r->s) swap(P->l,P->r);
	P->s = (P->r?P->r->s:0)+1; return P;
}
ph ins(ph p, T v) { return meld(p, new heap(v)); }
ph pop(ph p) { return meld(p->l,p->r); }
int N,M,src,des,K;
ph cand[MX];
vector<array<int,3>> adj[MX], radj[MX];
pi pre[MX];
ll dist[MX];
struct state {
	int vert; ph p; ll cost;
	bool operator<(const state& s) const { return cost > s.cost; }
};
int main() {
	setIO(); re(N,M,src,des,K);
	F0R(i,M) {
		int u,v,w; re(u,v,w);
		adj[u].pb({v,w,i}); radj[v].pb({u,w,i}); // vert, weight, label
	}
	priority_queue<state> ans;
	{
		F0R(i,N) dist[i] = INF, pre[i] = {-1,-1};
		priority_queue<T,vector<T>,greater<T>> pq;
		auto ad = [&](int a, ll b, pi ind) {
			if (dist[a] <= b) return;
			pre[a] = ind; pq.push({dist[a] = b,a});
		};
		ad(des,0,{-1,-1});
		vi seq;
		while (sz(pq)) {
			auto a = pq.top(); pq.pop(); 
			if (a.f > dist[a.s]) continue;
			seq.pb(a.s); trav(t,radj[a.s]) ad(t[0],a.f+t[1],{t[2],a.s}); // edge index, vert
		}
		trav(t,seq) {
			trav(u,adj[t]) if (u[2] != pre[t].f && dist[u[0]] != INF) {
				ll cost = dist[u[0]]+u[1]-dist[t];
				cand[t] = ins(cand[t],{cost,u[0]});
			}
			if (pre[t].f != -1) cand[t] = meld(cand[t],cand[pre[t].s]);
			if (t == src) {
				ps(dist[t]); K --;
				if (cand[t]) ans.push(state{t,cand[t],dist[t]+cand[t]->v.f});
			}
		}
	}
	F0R(i,K) {
		if (!sz(ans)) {
			ps(-1);
			continue;
		}
		auto a = ans.top(); ans.pop();
		int vert = a.vert;
		ps(a.cost);
		if (a.p->l) {
			ans.push(state{vert,a.p->l,a.cost+a.p->l->v.f-a.p->v.f});
		}
		if (a.p->r) {
			ans.push(state{vert,a.p->r,a.cost+a.p->r->v.f-a.p->v.f});
		}
		int V = a.p->v.s;
		if (cand[V]) ans.push(state{V,cand[V],a.cost+cand[V]->v.f});
	}
}

Pt LinesInter(Line a, Line b) {
    double abc = (a.b - a.a) ^ (b.a - a.a);
    double abd = (a.b - a.a) ^ (b.b - a.a);
    if (sign(abc - abd) == 0) return b.b;// no inter
    return (b.b * abc - b.a * abd) / (abc - abd);
}
vector <Pt> SegsInter(Line a, Line b) {
    if (btw(a.a, a.b, b.a)) return {b.a};
    if (btw(a.a, a.b, b.b)) return {b.b};
    if (btw(b.a, b.b, a.a)) return {a.a};
    if (btw(b.a, b.b, a.b)) return {a.b};
    if (ori(a.a, a.b, b.a) * ori(a.a, a.b, b.b) == -1 && ori(b.a, b.b, a.a) * ori(b.a, b.b, a.b) == -1) {
        return {LinesInter(a, b)};
    }
    return {};
}
double polyUnion(vector <vector <Pt>> poly) {
    int n = poly.size();
    double ans = 0;
    auto solve = [&](Pt a, Pt b, int cid) {
        vector <pair <Pt, int>> event;
        for (int i = 0; i < n; ++i) {
            int st = 0, sz = poly[i].size();
            while (st < sz && ori(poly[i][st], a, b) != 1) st++;
            if (st == sz) continue;
            for (int j = 0; j < sz; ++j) {
                Pt c = poly[i][(j + st) % sz], d = poly[i][(j + st + 1) % sz];
                if (sign((a - b) ^ (c - d)) != 0) {
                    int ok1 = ori(c, a, b) == 1;
                    int ok2 = ori(d, a, b) == 1;
                    if (ok1 ^ ok2) event.emplace_back(LinesInter({a, b}, {c, d}), ok1 ? 1 : -1);
                } else if (ori(c, a, b) == 0 && sign((a - b) * (c - d)) > 0 && i <= cid) {
                    event.emplace_back(c, -1);
                    event.emplace_back(d, 1);
                }
            }
        }
        sort(all(event), [&](pair <Pt, int> i, pair <Pt, int> j) {
            return ((a - i.first) * (a - b)) < ((a - j.first) * (a - b));
        });
        int now = 0;
        Pt lst = a;
        for (auto [x, y] : event) {
            if (btw(a, b, lst) && btw(a, b, x) && now == 0) ans += lst ^ x;
            now += y, lst = x;
        }
    };
    for (int i = 0; i < n; ++i) for (int j = 0; j < poly[i].size(); ++j) {
        Pt a = poly[i][j], b = poly[i][(j + 1) % int(poly[i].size())];
        solve(a, b, i);
    }
    return ans / 2;
}
// Minimum Steiner Tree, O(V 3^T + V^2 2^T)
struct SteinerTree { // 0-base
  static const int T = 10, N = 105, INF = 1e9;
  int n, dst[N][N], dp[1 << T][N], tdst[N];
  int vcost[N]; // the cost of vertexs
  void init(int _n) {
    n = _n;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) dst[i][j] = INF;
      dst[i][i] = vcost[i] = 0;
    }
  }
  void add_edge(int ui, int vi, int wi) {
    dst[ui][vi] = min(dst[ui][vi], wi);
  }
  void shortest_path() {
    for (int k = 0; k < n; ++k)
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
          dst[i][j] =
            min(dst[i][j], dst[i][k] + dst[k][j]);
  }
  int solve(const vector<int> &ter) {
    shortest_path();
    int t = SZ(ter);
    for (int i = 0; i < (1 << t); ++i)
      for (int j = 0; j < n; ++j) dp[i][j] = INF;
    for (int i = 0; i < n; ++i) dp[0][i] = vcost[i];
    for (int msk = 1; msk < (1 << t); ++msk) {
      if (!(msk & (msk - 1))) {
        int who = __lg(msk);
        for (int i = 0; i < n; ++i)
          dp[msk][i] =
            vcost[ter[who]] + dst[ter[who]][i];
      }
      for (int i = 0; i < n; ++i)
        for (int submsk = (msk - 1) & msk; submsk;
             submsk = (submsk - 1) & msk)
          dp[msk][i] = min(dp[msk][i],
            dp[submsk][i] + dp[msk ^ submsk][i] -
              vcost[i]);
      for (int i = 0; i < n; ++i) {
        tdst[i] = INF;
        for (int j = 0; j < n; ++j)
          tdst[i] =
            min(tdst[i], dp[msk][j] + dst[j][i]);
      }
      for (int i = 0; i < n; ++i) dp[msk][i] = tdst[i];
    }
    int ans = INF;
    for (int i = 0; i < n; ++i)
      ans = min(ans, dp[(1 << t) - 1][i]);
    return ans;
  }
};
