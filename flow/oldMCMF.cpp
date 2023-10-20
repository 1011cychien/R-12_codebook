template <typename Flow, typename Cost>
struct MinCostMaxFlow {
  static constexpr Flow flowINF = numeric_limits<Flow>::max() / 2;
  static constexpr Cost costINF = numeric_limits<Cost>::max() / 2;
  struct Edge {
    int to;
    Flow cap;
    Cost cost;
    Edge(int to, Flow cap, Cost cost) : to(to), cap(cap), cost(cost) {}
  };
  int n;
  vector<Edge> e;
  vector<vector<int>> g;
  vector<Cost> h, dis;
  vector<int> pre;
  MinCostMaxFlow(int n) : n(n), g(n) {}
  bool spfa(int s, int t) {
    dis.assign(n, costINF);
    pre.assign(n, -1);
    vector<int> q{s}, inq(n);
    dis[s] = 0;
    inq[s] = 1;
    for (int i = 0; i < int(q.size()); i++) {
      int u = q[i];
      inq[u] = 0;
      for (int j : g[u]) {
        auto [v, cap, cost] = e[j];
        if (Cost nd = dis[u] + cost; cap > 0 && nd < dis[v]) {
          dis[v] = nd;
          pre[v] = j;
          if (!inq[v]) {
            q.push_back(v);
            inq[v] = 1;
          }
        }
      }
    }
    return dis[t] != costINF;
  }
  bool dijkstra(int s, int t) {
    dis.assign(n, costINF);
    pre.assign(n, -1);
    priority_queue<pair<Cost, int>, vector<pair<Cost, int>>, greater<>> pq;
    dis[s] = 0;
    pq.emplace(0, s);
    while (!pq.empty()) {
      auto [d, u] = pq.top();
      pq.pop();
      if (dis[u] != d) continue;
      for (int i : g[u]) {
        auto [v, cap, cost] = e[i];
        if (Cost nd = d + h[u] - h[v] + cost; cap > 0 && dis[v] > nd) {
          dis[v] = nd;
          pre[v] = i;
          pq.emplace(dis[v], v);
        }
      }
    }
    return dis[t] != costINF;
  }
  void addEdge(int u, int v, Flow cap, Cost cost) {
    g[u].push_back(e.size());
    e.emplace_back(v, cap, cost);
    g[v].push_back(e.size());
    e.emplace_back(u, 0, -cost);
  }
  pair<Flow, Cost> maxFlow(int s, int t) {
    Flow flow = 0;
    Cost cost = 0;
    while (spfa(s, t)) {
      Flow aug = flowINF;
      for (int i = t; i != s; i = e[pre[i] ^ 1].to) {
        aug = min(aug, e[pre[i]].cap);
      }
      for (int i = t; i != s; i = e[pre[i] ^ 1].to) {
        e[pre[i]].cap -= aug;
        e[pre[i] ^ 1].cap += aug;
      }
      flow += aug;
      cost += aug * dis[t];
    }
    return make_pair(flow, cost);
  }
  pair<Flow, Cost> maxFlow2(int s, int t) {
    Flow flow = 0;
    Cost cost = 0;
    h.assign(n, 0);
    // Johnson's potential
    // Note that the graph must be DAG, and all edges must have u < v (S = 0, T = lst?)
    for (int i = 0; i < n; i++) {
      for (auto j : g[i]) {
        if (e[j].to > i) {
          h[e[j].to] = min(h[e[j].to], h[i] + e[j].cost);
        }
      }
    }
    while (dijkstra(s, t)) {
      for (int i = 0; i < n; ++i) {
        h[i] += dis[i];
      }
      Flow aug = flowINF;
      for (int i = t; i != s; i = e[pre[i] ^ 1].to) {
        aug = min(aug, e[pre[i]].cap);
      }
      for (int i = t; i != s; i = e[pre[i] ^ 1].to) {
        e[pre[i]].cap -= aug;
        e[pre[i] ^ 1].cap += aug;
      }
      flow += aug;
      cost += aug * h[t];
    }
    return make_pair(flow, cost);
  }
};
