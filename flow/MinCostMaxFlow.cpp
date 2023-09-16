template <class Flow, class Cost>
struct MinCostMaxFlow {
public:
    static constexpr Flow flowINF = numeric_limits<Flow>::max();
    static constexpr Cost costINF = numeric_limits<Cost>::max();
    MinCostMaxFlow() {}
    MinCostMaxFlow(int n) : n(n), g(n) {}
    int addEdge(int u, int v, Flow cap, Cost cost) {
        int m = int(pos.size());
        pos.push_back({u, int(g[u].size())});
        g[u].push_back({v, int(g[v].size()), cap, cost});
        g[v].push_back({u, int(g[u].size()) - 1, 0, -cost});
        return m;
    }
    struct edge {
        int u, v;
        Flow cap, flow;
        Cost cost;
    };
    edge getEdge(int i) {
        int m = int(pos.size());
        auto _e = g[pos[i].first][pos[i].second];
        auto _re = g[_e.v][_e.rev];
        return {pos[i].first, _e.v, _e.cap + _re.cap, _re.cap, _e.cost};
    }
    vector<edge> edges() {
        int m = int(pos.size());
        vector<edge> result(m);
        for (int i = 0; i < m; i++) { result[i] = getEdge(i); }
        return result;
    }
    pair<Flow, Cost> maxFlow(int s, int t, Flow flow_limit = flowINF) { return slope(s, t, flow_limit).back(); }
    vector<pair<Flow, Cost>> slope(int s, int t, Flow flow_limit = flowINF) {
        vector<Cost> dual(n, 0), dis(n);
        vector<int> pv(n), pe(n), vis(n);
        auto dualRef = [&]() {
            fill(dis.begin(), dis.end(), costINF);
            fill(pv.begin(), pv.end(), -1);
            fill(pe.begin(), pe.end(), -1);
            fill(vis.begin(), vis.end(), false);
            struct Q {
                Cost key;
                int u;
                bool operator<(Q o) const { return key > o.key; }
            };
            priority_queue<Q> h;
            dis[s] = 0;
            h.push({0, s});
            while (!h.empty()) {
                int u = h.top().u;
                h.pop();
                if (vis[u]) { continue; }
                vis[u] = true;
                if (u == t) { break; }
                for (int i = 0; i < int(g[u].size()); i++) {
                    auto e = g[u][i];
                    if (vis[e.v] || e.cap == 0) continue;
                    Cost cost = e.cost - dual[e.v] + dual[u];
                    if (dis[e.v] - dis[u] > cost) {
                        dis[e.v] = dis[u] + cost;
                        pv[e.v] = u;
                        pe[e.v] = i;
                        h.push({dis[e.v], e.v});
                    }
                }
            }
            if (!vis[t]) { return false; }
            for (int v = 0; v < n; v++) {
                if (!vis[v]) continue;
                dual[v] -= dis[t] - dis[v];
            }
            return true;
        };
        Flow flow = 0;
        Cost cost = 0, prevCost = -1;
        vector<pair<Flow, Cost>> result;
        result.push_back({flow, cost});
        while (flow < flow_limit) {
            if (!dualRef()) break;
            Flow c = flow_limit - flow;
            for (int v = t; v != s; v = pv[v]) {
                c = min(c, g[pv[v]][pe[v]].cap);
            }
            for (int v = t; v != s; v = pv[v]) {
                auto& e = g[pv[v]][pe[v]];
                e.cap -= c;
                g[v][e.rev].cap += c;
            }
            Cost d = -dual[s];
            flow += c;
            cost += c * d;
            if (prevCost == d) { result.pop_back(); }
            result.push_back({flow, cost});
            prevCost = cost;
        }
        return result;
    }
private:
    int n;
    struct _edge {
        int v, rev;
        Flow cap;
        Cost cost;
    };
    vector<pair<int, int>> pos;
    vector<vector<_edge>> g;
};