#include "BipartiteGraph.h"

#include <vector>
#include <cmath>
#include <algorithm>
#include <climits>
#include <iostream>
#include <unordered_map>
#include <cfloat>
#include <queue> 
#include <stack>
#include <set>
#include <utility>
#include <iomanip> 
#include <cassert>
const double eps = 1e-9;
const double INF = 1e20;

using namespace std;

void Alg1_bigraphCreate(const vector<pair<double,double>> &VehilcePosition, const vector<vector<int>> &kPath, const vector<pair<double,double>> &Rs, 
    const vector<pair<double,double>> &Rt, vector<vector<double>> &Adj){
    int X = VehilcePosition.size();

    auto  GetNewDist = [&](const int &u, const int &groupId){
        auto &[vx, vy] = VehilcePosition[u];
        auto &vec = kPath[groupId];

        double MinDist = 100000000;
        for(const auto &sid: vec){
            auto &[sx, sy] = Rs[sid];
            auto &[tx, ty] = Rt[sid];
            double dis_us = sqrt((vx - sx) * (vx - sx) + (vy - sy)*(vy - sy));
            double dis_st = sqrt((sx - tx) * (sx - tx) + (sy - ty)*(sy - ty));
            MinDist = min(MinDist, dis_us + dis_st);
        }
        return MinDist;
    };
    Adj.resize(X);
    for(int i = 0; i < X; ++ i){ // vehicle
        Adj[i].resize(X);
        for(int j = 0; j < X; ++ j){ // groupid
            double dist = GetNewDist(i, j);
            Adj[i][j] = dist;
        }
    }
}

void Alg2_bigraphCreate(const vector<pair<double,double>> &VehilcePosition, const vector<vector<int>> &kPath, const vector<pair<double,double>> &Rs, 
    const vector<pair<double,double>> &Rt, vector<vector<double>> &Adj, const int &gap){
    int X = VehilcePosition.size();

    auto  GetNewDist = [&](const int &u, const int &groupId){
        auto &[vx, vy] = VehilcePosition[u];
        auto &vec = kPath[groupId];

        double MinDist = 100000000;
        for(const auto &sid: vec){
            // the number of t will bigger than gap
            if(sid > gap) continue;
            auto &[sx, sy] = Rs[sid];
            auto &[tx, ty] = Rt[sid];
            double dis_us = sqrt((vx - sx) * (vx - sx) + (vy - sy)*(vy - sy));
            double dis_st = sqrt((sx - tx) * (sx - tx) + (sy - ty)*(sy - ty));
            MinDist = min(MinDist, dis_us + dis_st);
        }
        return MinDist;
    };
    Adj.resize(X);
    for(int i = 0; i < X; ++ i){ // vehicle
        Adj[i].resize(X);
        for(int j = 0; j < X; ++ j){ // groupid
            double dist = GetNewDist(i, j);
            Adj[i][j] = dist;
        }
    }
}

bool KM::dfs(int u){
    va[u] = true;
    for (int v = 0; v < n; v++) {
        if (!vb[v]) {
            double diff = la[u] + lb[v] - w[u][v];
            if (fabs(diff) < eps) {  // 
                vb[v] = true;
                if (match[v] == -1 || dfs(match[v])) {
                    match[v] = u;
                    return true;
                }
            } else if (slack[v] - diff > eps) {  // slack[v] > diff
                slack[v] = diff;
            }
        }
    }
    return false;
}

double KM::solve(){

    for(int i = 0; i < n; ++ i){
        for(int j = 0; j < n; ++ j){
            la[i] = max(la[i], w[i][j]);
        }
    }
    for (int i = 0; i < n; i++) {
        slack.assign(n, INF);
        
        while (true) {
            va.assign(n, false);
            vb.assign(n, false);
            
            if (dfs(i)) break;
            
            // compute delta
            double delta = INF;
            for (int j = 0; j < n; j++) {
                if (!vb[j] && slack[j] < delta) {
                    delta = slack[j];
                }
            }
            // adjust
            for (int j = 0; j < n; j++) {
                if (va[j]) la[j] -= delta;
                if (vb[j]) lb[j] += delta;
                else slack[j] -= delta;
            }
        }
    }

    double ans = 0;
    for (int j = 0; j < n; j++) {
        ans += w[match[j]][j];
    }
    return -ans;  
    return 0;

}

void KM::init(int n){
    this->n = n;
    la.resize(n, -INF);
    lb.resize(n, 0);
    match.resize(n, -1);
    w.resize(n, vector<double>(n));
}

void VehicleminCostMatching(const vector<vector<double>> &Adj, vector<int> &Result) {

    int n = Adj.size();
    Result.resize(n);
    KM km;
    km.init(n);
    for(int i = 0; i < n; ++ i){
        for(int j = 0; j < n; ++ j){
            km.w[i][j] = -Adj[i][j];        
        }
    }

    double ans = km.solve();
    printf("VehicleminCost is : %.6lf\n", ans);
    for(int j = 0; j < n; ++ j){
        cout << "Vehicle " << km.match[j] << " choose group is : " << j << '\n';
        Result[km.match[j]] = j;
    }
}

void GroupminCostMatching(const vector<vector<double>> &Adj, vector<int> &Result) {

    int n = Adj.size();
    Result.resize(n);

    KM km;
    km.init(n);
    for(int i = 0; i < n; ++ i){
        for(int j = 0; j < n; ++ j){
            km.w[i][j] = -Adj[i][j];        
        }
    }

    double ans = km.solve();
    // printf(" %.6lf\n", ans);
    cout << " GroupminCost is :" << fixed << setprecision(6) << ans << endl;
    for(int j = 0; j < n; ++ j){
        cout << "Group " << km.match[j] << " choose group is : " << j << '\n';
        Result[km.match[j]] = j;
    }
}

static void dfsEuler(int u, std::unordered_map<int, std::vector<int>>& adj, std::vector<int>& walk) {
    std::stack<int> stk;
    stk.push(u);
    while (!stk.empty()) {
        int v = stk.top();
        if (adj[v].empty()) {
            walk.push_back(v);
            stk.pop();
        } else {
            int u2 = adj[v].back();
            adj[v].pop_back();
            stk.push(u2);
        }
    }
    std::reverse(walk.begin(), walk.end());
}

void Alg1_computeResult(const vector<pair<double,double>> &VehilcePosition, const vector<vector<int>> &kPath, const vector<pair<double,double>> &Rs, 
    const vector<pair<double,double>> &Rt, const vector<int> &Result){
    
    int m = VehilcePosition.size();
    int k = kPath[0].size();

    auto GetEdgeVal = [&](pair<double,double> xa, pair<double,double> xb){
        double dist_x = sqrt((xa.first - xb.first) * (xa.first - xb.first) + (xa.second - xb.second) * (xa.second - xb.second));
        return dist_x;
    };

    double sum = 0;

    for(int vehicleId = 0; vehicleId < m; ++ vehicleId){ 
        int groupiId = Result[vehicleId]; 
        const auto &path = kPath[groupiId];
        const auto &vpos = VehilcePosition[vehicleId];
        cout << "**************************\n";
        cout << "Vehicle Id is : " << vehicleId << endl;
        cout << "Group Id is :" << groupiId << endl;
        cout << "The path is : ";

        for(const auto &it: path){
            cout << it << ' ';
        }
        cout << endl;
        cout << "**************************\n";


        vector<int> EulerPath;
        double MinCost = INF;
        // u define as  the point 2*k + vehicleId
        int u = 2*m*k+vehicleId;

        const int &lefts = path[0];
        const int &rights = path[k-1];
        const int &leftt = path[0] + m*k;
        const int &rightt = path[k-1] + m*k;

        // path[i] = {si, ti}
        // si --> path[i], ti --> path[i] + m*k;
        for(int i = 0; i < k; ++ i){ // select path[i] as the start point
            // addedge(u, si)
            std::unordered_map<int, std::vector<int>> adj_eluer;
            
            for(int j = 1; j < k; ++ j){
                int su = path[j-1];
                int sv = path[j];
                adj_eluer[su].push_back(sv);
                adj_eluer[sv].push_back(su);
            }
           
            const int &si = path[i];

            adj_eluer[u].push_back(si);
            adj_eluer[si].push_back(si + m*k);

            const double t_left = GetEdgeVal(Rt[lefts], Rt[si]);
            const double t_right = GetEdgeVal(Rt[rights], Rt[si]);
            
            if(t_left < t_right){
                if(si+m*k != leftt){
                    adj_eluer[si + m*k].push_back(leftt);  
                }          
                for(int j = 1; j < k; ++ j){
                    int tu = path[j-1] + m*k;
                    int tv = path[j] + m*k;
                    adj_eluer[tu].push_back(tv);
                }
            }else{
                if(rightt != si + m*k){
                    adj_eluer[si + m*k].push_back(rightt);  
                }
                for(int j = 1; j < k; ++ j){
                    int tu = path[j-1] + m*k;
                    int tv = path[j] + m*k;
                    adj_eluer[tv].push_back(tu);
                }
            }
            // test
            // g.printGraph();
            if (1) {
                cout << "\nEulerian Path Exists!" << endl;             
                vector<int> temp;
                dfsEuler(u, adj_eluer, temp);
                cout << "Path: ";
                int last = -1;
                for (int node : temp) {
                    cout << node << " ";
                }
                cout << endl;

                set<int> se;
                vector<int> Epath;
                double tot = 0;
                cout << "Remove the same nodes:" << endl;
                pair<double,double> last_pos;
                bool ff = false;
                for(const auto &v: temp){
                    if(se.count(v)) continue;
                    se.insert(v);
                    Epath.push_back(v);
                    int op, idx;
                    pair<double,double> now;
                    if(v < m * k){
                        cout << "S" << v << "->";
                        idx = v;
                        now = Rs[idx];
                    }else if(v < 2 * m *k){
                        idx = v - m*k;
                        cout << "T" << idx << "->";
                        now = Rt[idx];
                    }else{
                        idx = v - 2 *m*k;
                        cout << "V" << idx << "->";
                        now = vpos;
                    }
                    if(ff) {
                        tot += GetEdgeVal(last_pos, now);
                    }
                    last_pos = now;
                    ff = true;
                }
                cout << endl;
                if(tot < MinCost){
                    MinCost = tot;
                    EulerPath = Epath;
                }
            } else {
                cout << "\nNo Eulerian Path Exists!" << endl;
            }
        }
        
        sum += MinCost;
        cout << "----------------------" << endl;
        cout << "Vehicle number is :" << vehicleId << endl;
        cout << "minumum Cost is : " << fixed << setprecision(6) << MinCost << endl;
        cout << "The minumum cost Path is: \n";
        for(const auto &u: EulerPath){
            if(u < m * k){
                cout << "S" << u << "->";
            }else if(u < 2 * m *k){
                    int idx = u - m*k;
                    cout << "T" << idx << "->";
            }else{
                int idx = u - 2 *m*k;
                cout << "V" << idx << "->";
            }
        }
        cout << '\n';
        cout << "----------------------" << endl;
    }
    cout << "TOT COST is : " << fixed << setprecision(6)  << sum << '\n';
}

void Alg2_computeResult(const vector<pair<double,double>> &VehilcePosition, const vector<vector<int>> &kPath, const vector<pair<double,double>> &Rs, 
    const vector<pair<double,double>> &Rt, const vector<int> &Result, const int &gap)
    {

    
    int m = VehilcePosition.size();
    // kpath is a complete path, consist set s and set t
    int k = kPath[0].size() / 2;
    double tot_cost = 0;
    auto GetEdgeVal = [&](pair<double,double> xa, pair<double,double> xb){
        double dist_x = sqrt((xa.first - xb.first) * (xa.first - xb.first) + (xa.second - xb.second) * (xa.second - xb.second));
        return dist_x;
    };

    for(int vehicleId = 0; vehicleId < m; ++ vehicleId){ 
        int groupiId = Result[vehicleId]; 
        const auto &path = kPath[groupiId];
        const auto &vpos = VehilcePosition[vehicleId];
        cout << "**************************\n";
        cout << "Vehicle Id is : " << vehicleId << endl;
        cout << "Group Id is :" << groupiId << endl;
        cout << "The path is : ";

        for(const auto &it: path){
            cout << it << ' ';
        }
        cout << endl;
        cout << "**************************\n";
        
        // Graph eulerGraph(2*m*k + m);

        vector<int> EulerPath;
        double MinCost = INF;
        // u define as  the point 2*k + vehicleId
        int u = 2*m*k+vehicleId;
        vector<int> set_s, set_t;
        set<int> ss;
        for(int i = 0; i < 2*k; ++ i){
            if(path[i] < m*k) {
                set_s.push_back(path[i]);
                ss.insert(path[i]);
            }
            else {
                assert(ss.count(path[i] - m*k));
                set_t.push_back(path[i]);
            }
        }


        const int &lefts = set_s[0];
        const int &rights = set_s[k-1];
        const int &leftt = set_t[0];
        const int &rightt = set_t[k-1];

        for(int i = 0; i < k; ++ i){
            std::unordered_map<int, std::vector<int>> g;
            for(int j = 1; j < k; ++ j){
                g[set_s[j-1]].push_back(set_s[j]);
                g[set_s[j]].push_back(set_s[j-1]);
            }
            // Graph g = eulerGraph;

            const auto &s_now = set_s[i];
            const auto &t_now = set_t[i];

            g[s_now].push_back(t_now);

            g[u].push_back(s_now);

            const double &dist_lt = GetEdgeVal(Rt[leftt - gap], Rt[t_now - gap]);
            const double &dist_rt = GetEdgeVal(Rt[t_now - gap], Rt[rightt - gap]);
            if(dist_lt - dist_rt <= eps){
                if(leftt != t_now){
                    g[t_now].push_back(leftt);  
                }
                for(int j = 1; j < k; ++ j){
                    g[set_t[j-1]].push_back(set_t[j]);
                }
            }else{
                if(rightt != t_now){
                    g[t_now].push_back(rightt);  
                }
                for(int j = 1; j < k; ++ j){
                    g[set_t[j]].push_back(set_t[j-1]);
                }
            }
            // graph construct complete

            if (1) {
                cout << "\nEulerian Path Exists!" << endl;    
                cout << "u to s_i: " <<  vehicleId << " -> " << s_now << endl;        
                // auto temp = g.findEulerPath(u);
                vector<int> temp;
                dfsEuler(u, g, temp);
                cout << "Path: ";
                int last = -1;
                for (int node : temp) {
                    cout << node << " ";
                }
                cout << endl;

                set<int> se;
                vector<int> Epath;
                double tot = 0;
                cout << "Remove the same nodes:" << endl;
                pair<double,double> last_pos;
                bool ff = false;
                for(const auto &v: temp){
                    if(se.count(v)) continue;
                    se.insert(v);
                    Epath.push_back(v);
                    int op, idx;
                    pair<double,double> now;
                    if(v < gap){
                        cout << "S" << v << "->";
                        idx = v;
                        now = Rs[idx];
                    }else if(v < 2 * m * k){
                        idx = v - m*k;
                        cout << "T" << idx << "->";
                        now = Rt[idx];
                    }else{
                        idx = v - 2 *m*k;
                        cout << "V" << idx << "->";
                        now = vpos;
                    }
                    if(ff) {
                        tot += GetEdgeVal(last_pos, now);
                    }
                    last_pos = now;
                    ff = true;
                }
                cout << endl;
                if(tot - MinCost < 1e-3){
                    MinCost = tot;
                    EulerPath = Epath;
                }
                cout << "This case's cost is: " << fixed << setprecision(6) <<  tot << endl;

            } else {
                cout << "\nNo Eulerian Path Exists!" << endl;
            }
        }


        cout << "----------------------" << endl;
        cout << "Vehicle number is :" << vehicleId << endl;
        cout << "minumum Cost is : " << fixed << setprecision(6) << MinCost << endl;
        cout << "The minumum cost Path is: \n";
        for(const auto &u: EulerPath){
            if(u < gap)
                cout << "S" << u << "->" ;
            else if(u < 2 * m * k)
                cout << "T" << u - gap << "->" ;
            else 
                cout << "V" << u - 2 * m * k << "->" ;
        }
        cout << '\n';
        cout << "----------------------" << endl;

        tot_cost += MinCost;
    }
    
    cout << "TOT COST is: " << tot_cost << fixed << setprecision(6) << endl;
}