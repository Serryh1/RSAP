#include "SF-Solution.h"
#include "SF-Functions.h"
#include "GW.h"

#define __STDC_LIMIT_MACROS
#include <iostream>
#include <time.h>
#include <sys/time.h>
#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <set>
#include <limits>
#include <queue>
#include <cassert>
#include <algorithm>
#include <functional>
#include <numeric>
#include <queue>
#include <set>
#include <tuple>
#include <unordered_set>
#include <vector>
#include <map>
#include <chrono>
using namespace std;

class edges
{
public:
    double edgeWeight;
    pair<int,int> edgeCoordinates;
};

bool operator< (const edges& lhs, const edges& rhs)
{
    return lhs.edgeWeight < rhs.edgeWeight;
}

bool operator> (const edges& lhs, const edges& rhs)
{
    return lhs.edgeWeight > rhs.edgeWeight;
}

bool operator== (const edges& lhs, const edges& rhs)
{
    return (lhs.edgeCoordinates.first == rhs.edgeCoordinates.first)||(lhs.edgeCoordinates.second == rhs.edgeCoordinates.first);
}

vector<int> father;
vector<int> _size;
vector<int> Rank;
vector<int> isactive;
int union_find(int x) {
    if (father[x] != x) father[x] = union_find(father[x]);
    return father[x];
}

void unite(int x, int y, int K) {
    int px = union_find(x);
    int py = union_find(y);
    if (px == py) return;

    if (Rank[px] < Rank[py]) {
        father[px] = py;
        _size[py] += _size[px];
        if(_size[py] % K == 0) {
            isactive[py] = 0;
        }else{
            isactive[py] = 1;
        }
    } else if (Rank[px] > Rank[py]) {
        father[py] = px;
        _size[px] += _size[py];
        if(_size[px] % K == 0) {
            isactive[px] = 0;
        }else{
            isactive[px] = 1;
        }
    } else {
        father[py] = px;
        _size[px] += _size[py];
        if(_size[px] % K == 0) {
            isactive[px] = 0;
        }else{
            isactive[px] = 1;
        }
        Rank[px]++;
    }
}


void unite1(int x, int y, int K) {
    int px = union_find(x);
    int py = union_find(y);
    if (px == py) return;

    if (Rank[px] < Rank[py]) {
        father[px] = py;
        _size[py] += _size[px];
    } else if (Rank[px] > Rank[py]) {
        father[py] = px;
        _size[px] += _size[py];
    } else {
        father[py] = px;
        _size[px] += _size[py];
        Rank[px]++;
    }
}

struct Component {
    std::vector<int> nodes;
    std::vector<std::pair<int, int>> edges;
    std::vector<int> degree;
    double dual = 0.0;
    bool frozen = false;

    Component(int n, int id) : degree(n, 0) { nodes.push_back(id); }
};

struct Edge{

    double key;
    int u, v;
    int label;
    Edge(double k, int u, int v, int label) : key(k), u(u), v(v), label(label) {}

    bool operator>(const Edge& other) const {
        return key > other.key;
    }
};

struct newedge{
    // Component (id, node)
    pair<int,int> Componentu, Componentv;
    int label;
    
    bool friend operator<(newedge x, newedge y){
        return x.Componentv.first < y.Componentv.first;
    };
};

void GW_ALG1(int &nVertex, int &nEdges, int &nTerminals, int *&terminals, vector<vector<double > > & adjMatrix, 
    vector<vector<int > > & pairedTerminals, vector<vector<int > > & E , const int &K, vector<vector<int>> &kPath)
{   

    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();

    steinerForest solutionGraph(nVertex); // Create the solution graph
    solutionGraph.Path.clear(); // path
    solutionGraph.visit.resize(nVertex, false);

    // for clean process
    for(const auto &it: pairedTerminals){
        solutionGraph.Terminals[it[0]] = it[0];
    }

    cout << "K is" << K << endl;
    vector<double> d(nVertex);

    // union initial
    father.resize(nVertex, 0);
    _size.resize(nVertex, 1);
    Rank.resize(nVertex, 0);
    isactive.resize(nVertex, 1);

    int n = nVertex;
    for (int i = 0; i < n; ++i) father[i] = i;

    // store the edges
    vector<pair<int,int>> solutionsEdges;
    priority_queue<Edge, vector<Edge>, greater<Edge>> pq;

    vector<vector<newedge>> new_Edges(n);
    // for each components build an edge set <v, label>
    for(int i = 0; i < n; ++ i){
        for(int j = i + 1; j < n; ++ j){
            pq.push(Edge(adjMatrix[i][j] / 2, i, j, 0));
            pair<int,int> C_u = {i, i};
            pair<int,int> C_v = {j, j};
            new_Edges[i].push_back({C_u, C_v, 0});
            new_Edges[j].push_back({C_v, C_u, 0});
        }
        sort(new_Edges[i].begin(), new_Edges[i].end());
    }

    bool flag = false;
    // check there is no active nodes
    tuple<int,int,int> tp;

    // restore the deleted edges
    map<tuple<int,int,int>,int> LazyPQ;
    double last_key = 0;
    int idx = 0;
    // the label idx
    while(!pq.empty() && !flag){
        idx ++;
        const auto &[key, u, v, label] = pq.top();
        pq.pop();
        // if this edge has been deleted, then continue
        if(LazyPQ.count(make_tuple(u, v, label))) continue;
        int Cp = union_find(u), Cq = union_find(v);
        // if the two components are inactive or they belong to the same component, then continue
        if(Cp == Cq || (!isactive[Cp] && !isactive[Cq])) continue;

        solutionGraph.addEdge(u, v, adjMatrix[u][v]);
        solutionsEdges.push_back({u, v});

        double epsilon  = key - last_key;
        last_key = key;
        // undate all active nodes
        for(int i = 0; i < n; ++ i){
            int p_i = union_find(i);
            if(isactive[p_i]){
                d[i] += epsilon;
            }
        }
        // union components
        if((_size[Cp] + _size[Cq]) % K == 0){
            isactive[Cp] = 0;
            isactive[Cq] = 0;
        }else{
            isactive[Cp] = 1;
            isactive[Cq] = 1;
        }

        // delete edges

        int size_cp = new_Edges[Cp].size();
        int size_cq = new_Edges[Cq].size();
        int j = 0;

        vector<newedge> temp;

        for(int i = 0; i < size_cp; ++ i){
            const auto &[ComponentP, ComponentPR, labelp] = new_Edges[Cp][i];
            while(j < size_cq && new_Edges[Cq][j].Componentv.first < ComponentPR.first){
                j ++;
            }
            if(j == size_cq) continue;
            auto [ComponentQ, ComponentQR, labelq] = new_Edges[Cq][j];
            // beglongs to the same 

            int pu = ComponentP.second;
            int pr = ComponentPR.second;
            int qu = ComponentQ.second;
            int qr = ComponentQR.second;
            if(ComponentPR.first != ComponentQR.first) continue;
            // delete two edges 
            LazyPQ[make_tuple(pu, pr, labelp)] = 1;
            LazyPQ[make_tuple(qu, qr, labelq)] = 1;
            double key_p = key;
            double key_q = key;

            int Cr = ComponentPR.first;
            
            int activepr = isactive[Cp] + isactive[Cr];
            int activeqr = isactive[Cq] + isactive[Cr];
            
            if(!activepr && !activeqr) {
                double p_1 = adjMatrix[pu][pr] - d[pu] - d[pr];
                double p_2 = adjMatrix[qu][qr] - d[qu] - d[qr];

                if(p_1 - p_2 <= 1e-3){
                    temp.push_back({ComponentP, ComponentPR, idx});
                    pq.push({1e9, pu, pr, idx});
                }else{
                    temp.push_back({ComponentQ, ComponentQR, idx});
                    pq.push({1e9, qu, qr, idx});
                }
                continue;
            }

            
            key_p += (adjMatrix[pu][pr] - d[pu] - d[pr]) / activepr;
            key_q += (adjMatrix[qu][qr] - d[qu] - d[qr]) / activeqr;
            
            // add back new min key edges
            if(key_p - key_q <= 1e-3){
                temp.push_back({ComponentP, ComponentPR, idx});
                pq.push({key_p, pu, pr, idx});
            }else{  
                temp.push_back({ComponentQ, ComponentQR, idx});
                pq.push({key_q, qu, qr, idx});
            }
        }
        new_Edges[Cp].clear();
        new_Edges[Cq].clear();

        unite(Cp, Cq, K);
        int C_new = union_find(Cp);
        sort(temp.begin(), temp.end());
        
        for(auto &[ComP, ComQ, label]: temp){
            ComP.first = C_new;
            new_Edges[C_new].push_back({ComP, ComQ, label});
        }
    }   

    cout<<"Out of while loop \n";

    // BFS Clean process
    for(int i = 0; i < solutionsEdges.size(); ++ i){
        auto [u, v] = solutionsEdges[i];
        if(solutionGraph.isvalidcut(u, v, nVertex / 2, K)){
            solutionGraph.deleteEdge(u, v);
        }
    }

    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "GW running time is: " << duration.count() << " s" << std::endl;


    vector<bool> visited(nVertex); 
    int counter = 0;
    // cout << K << '\n';
    // group split process
    for(int i = 0; i < nVertex; ++ i){
        if(!visited[i]){
            cout << "Counter is:" << counter << '\n';
            cout << "The path is: ";
            solutionGraph.Path.clear();
            solutionGraph.dfs(i);
            cout << '\n';
            counter ++;
            int cnt = 0;
            vector<int> temp;
            for(const int &u: solutionGraph.Path){
                cout << u << ' ';
                cnt ++;
                visited[u] = true;
                temp.push_back(u);
                if(cnt == K){
                    kPath.push_back(temp);
                    temp.clear();
                    cnt = 0;
                    cout << '\n';
                }
            }
            cout << '\n';
        }
    }
};


void sequenceSplit(const vector<int> &s, const vector<int> &t, vector<vector<int>> &kPath,const int &k, const int &gap){

    // this k is K / 2, only obtain a group of s or the t
    queue<vector<int>> seq_s, seq_t;
    seq_s.push(s);
    seq_t.push(t);

    while(!seq_s.empty()){
        auto vec_s = seq_s.front();
        seq_s.pop();

        int len = vec_s.size();

        // if len < 4 * k, we split the sequence using by brute force
        if(len < 4 * k){
            assert(len % k == 0);
            auto vec_t = seq_t.front();
            seq_t.pop();

            vector<vector<int>> temp_s;
            vector<int> vec;
            // split s
            unordered_map<int,int> label; // each si belongs to label_i
            for(int i = 1; i <= len; ++ i){
                label[vec_s[i-1]] = (i + k - 1) / k;
                vec.push_back(vec_s[i-1]);
                if(i % k == 0){
                    temp_s.push_back(vec);
                    vec.clear();
                }
            }
            int group_cnt = (len + k - 1) / k;
            vector<vector<int>> temp_t;
            temp_t.resize(group_cnt);

            for(const auto &tt: vec_t){
                assert(label.count(tt - gap));
                int lb = label[tt - gap];
                temp_t[lb-1].push_back(tt);
            }
            for(int i = 0; i < group_cnt; ++ i){
                vec.clear();
                const auto &vs = temp_s[i];
                const auto &vt = temp_t[i];
                for(const auto &v: vs){
                    vec.push_back(v);
                }
                for(const auto &v: vt){
                    vec.push_back(v);
                }
                kPath.push_back(vec);
            }
            continue;
        }

        // The first split process
        int left = (len + 2*k - 1) / (2*k) * k;
        int right = len - left;
        vector<int> ls, rs;
        unordered_map<int,int> label1;
        for(int i = 0; i < len; ++ i){
            if(i < left) {
                ls.push_back(vec_s[i]);
                label1[vec_s[i]] = 1;
            }else{
                rs.push_back(vec_s[i]);
                label1[vec_s[i]] = 2;
            }
        }

        const auto vec_t = seq_t.front();
        seq_t.pop();
        vector<int> lt, rt;
        for(const auto &v: vec_t){
            assert(label1.count(v - gap));
            if(label1[v - gap] == 1) {
                lt.push_back(v);
            }else{
                rt.push_back(v);
            }
        }
        // The second split 
        // split lt
        int ltsize = lt.size();
        int lpos = (ltsize + 2*k - 1) / (2*k) * k;
        vector<vector<int>> temp_t(4), temp_s(4);
        unordered_map<int,int> label2;
        for(int i = 0; i < ltsize; ++ i){
            int v = lt[i];
            if(i < lpos){
                label2[v - gap] = 0;
                temp_t[0].push_back(v);
            }else{
                label2[v - gap] = 1;
                temp_t[1].push_back(v);
            }
        }
        // split rt
        int rtsize = rt.size();
        int rpos = (rtsize + 2*k - 1) / (2*k) * k;
        for(int i = 0; i < rtsize; ++ i){
            int v = rt[i];
            if(i < rpos){
                label2[v - gap] = 2;
                temp_t[2].push_back(v);
            }else{
                label2[v - gap] = 3;
                temp_t[3].push_back(v);
            }
        }

        // split ls
        for(const auto &v: ls){
            assert(label2.count(v));
            int c = label2[v];
            temp_s[c].push_back(v);
        }
        // split rs
        for(const auto &v: rs){
            assert(label2.count(v));
            int c = label2[v];
            temp_s[c].push_back(v);
        }

        unordered_map<int,int> mp;
        for(int i = 0; i < 4; ++ i){
            const auto &vs = temp_s[i];
            const auto &vt = temp_t[i];
            
            for(const auto &v: vs){
                mp[v] = 1;
            }
            for(const auto &v: vt){
                assert(mp[v-gap] == 1);
            }
            seq_s.push(vs);
            seq_t.push(vt);
        }
    }

}


void GW_ALG2(int &nVertex, int &nEdges, int &nTerminals, int *&terminals, vector<vector<double > > & adjMatrix, 
    vector<vector<int > > & pairedTerminals, vector<vector<int > > & E , const int &K, vector<vector<int>> &kPath)
{

    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();

    cout << "nVertex is: " << nVertex << endl;
    cout << "K is " << K << endl;
    // 2*n vertices
    steinerForest solutionGraph(nVertex); // Create the solution graph
    solutionGraph.Path.clear(); // path
    solutionGraph.visit.resize(nVertex, false);

    for(const auto &it: pairedTerminals){
        solutionGraph.Terminals[it[0]] = it[1];
        solutionGraph.Terminals[it[1]] = it[0];
    }

    // growth 
    vector<double> d(nVertex);

    father.resize(nVertex, 0);
    _size.resize(nVertex, 1);
    Rank.resize(nVertex, 0);
    isactive.resize(nVertex, 1);

    vector<pair<int,int>> solutionsEdges;

    int n = nVertex;
    for (int i = 0; i < n; ++i) father[i] = i;

    priority_queue<Edge, vector<Edge>, greater<Edge>> pq;

    vector<vector<newedge>> new_Edges(n);
    // for each components build an edge set <v, label>
    for(int i = 0; i < n; ++ i){
        for(int j = i + 1; j < n; ++ j){
            pq.push(Edge(adjMatrix[i][j] / 2, i, j, 0));
            pair<int,int> C_u = {i, i};
            pair<int,int> C_v = {j, j};
            new_Edges[i].push_back({C_u, C_v, 0});
            new_Edges[j].push_back({C_v, C_u, 0});
        }
        sort(new_Edges[i].begin(), new_Edges[i].end());
    }

    bool flag = false;
    // check there is no active nodes

    tuple<int,int,int> tp;

    map<tuple<int,int,int>,int> LazyPQ;
    // restore the deleted edges
    double last_key = 0;
    int idx = 0;
    while(!pq.empty() && !flag){
        idx ++;
        auto [key, u, v, label] = pq.top();
        
        pq.pop();
        // if this edge has been deleted, then continue
        if(LazyPQ.count(make_tuple(u, v, label))) continue;
        int Cp = union_find(u), Cq = union_find(v);
        // if the two components are inactive
        if(Cp == Cq || (!isactive[Cp] && !isactive[Cq])) continue;

        solutionGraph.addEdge(u, v, adjMatrix[u][v]);
        solutionsEdges.push_back(make_pair(u, v));

        double epsilon  = key - last_key;
        last_key = key;
        // undate all active nodes
        for(int i = 0; i < n; ++ i){
            int p_i = union_find(i);
            if(isactive[p_i]){
                d[i] += epsilon;
            }
        }
        // cout << "u v is: " << u << ' ' << v << endl;
        // cout << "epsilon is :" << epsilon << endl;
        // union components

        int x;
        auto check = [&](int com1, int com2){
            set<int> vec;
            unordered_map<int,int> mp;
            for(const auto &[comp, comr, label]: new_Edges[com1]){
                vec.insert(comp.second);
                mp[comp.second] = 1;
            }
            for(const auto &[comp, comr, label]: new_Edges[com2]){
                vec.insert(comp.second);
                mp[comp.second] = 1;
            }
            
            for(const auto it: vec){
                // 
                if(mp[it + nVertex / 2] == 1 || mp[it - nVertex / 2] == 1){
                    continue;
                }else{
                    return false;
                }
            }
            return true;
        };

        if((_size[Cp] + _size[Cq]) % K == 0 && check(Cp, Cq)){
            isactive[Cp] = 0;
            isactive[Cq] = 0;
        }else{
            isactive[Cp] = 1;
            isactive[Cq] = 1;
        }
        // delete edges

        int size_cp = new_Edges[Cp].size();
        int size_cq = new_Edges[Cq].size();
        int j = 0;

        vector<newedge> temp;

        for(int i = 0; i < size_cp; ++ i){
            auto [ComponentP, ComponentPR, labelp] = new_Edges[Cp][i];
            while(j < size_cq && new_Edges[Cq][j].Componentv.first < ComponentPR.first){
                j ++;
            }
            if(j == size_cq) continue;
            auto [ComponentQ, ComponentQR, labelq] = new_Edges[Cq][j];
            // beglongs to the same 

            int pu = ComponentP.second;
            int pr = ComponentPR.second;
            int qu = ComponentQ.second;
            int qr = ComponentQR.second;
            if(ComponentPR.first != ComponentQR.first) continue;
            // delete two edges 
            LazyPQ[make_tuple(pu, pr, labelp)] = 1;
            LazyPQ[make_tuple(qu, qr, labelq)] = 1;
            double key_p = key;
            double key_q = key;
            int Cr = ComponentPR.first;
            
            int activepr = isactive[Cp] + isactive[Cr];
            int activeqr = isactive[Cq] + isactive[Cr];
            
            if(!activepr && !activeqr) {
                double p_1 = adjMatrix[pu][pr] - d[pu] - d[pr];
                double p_2 = adjMatrix[qu][qr] - d[qu] - d[qr];

                if(p_1 - p_2 <= 1e-3){
                    temp.push_back({ComponentP, ComponentPR, idx});
                    pq.push({1e9, pu, pr, idx});
                }else{
                    temp.push_back({ComponentQ, ComponentQR, idx});
                    pq.push({1e9, qu, qr, idx});
                }
                continue;
            }
            key_p += (adjMatrix[pu][pr] - d[pu] - d[pr]) / activepr;
            key_q += (adjMatrix[qu][qr] - d[qu] - d[qr]) / activeqr;
            

            // add back new min key edges
            if(key_p - key_q <= 1e-3){
                temp.push_back({ComponentP, ComponentPR, idx});
                pq.push({key_p, pu, pr, idx});
            }else{  
                temp.push_back({ComponentQ, ComponentQR, idx});
                pq.push({key_q, qu, qr, idx});
            }
        }
        new_Edges[Cp].clear();
        new_Edges[Cq].clear();

        unite1(Cp, Cq, K);
        int C_new = union_find(Cp);
        sort(temp.begin(), temp.end());
        new_Edges[C_new].clear();
        for(auto &[ComP, ComQ, label]: temp){
            ComP.first = C_new;
            new_Edges[C_new].push_back({ComP, ComQ, label});
        }
    }   

    // BFS Clean process
    for(int i = 0; i < solutionsEdges.size(); ++ i){
        auto [u, v] = solutionsEdges[i];
        if(solutionGraph.isvalidcut(u, v, nVertex / 2, K)){
            solutionGraph.deleteEdge(u, v);
        }
    }
    
    cout<<"Out of while loop \n";
        
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "GW running time is: " << duration.count() << " s" << std::endl;


    vector<bool> visited(nVertex); 
    int counter = 0;
    cout << K << '\n';
    for(int i = 0; i < nVertex; ++ i){
        if(!visited[i]){
            cout << "Counter :" << counter  << '\n';
            cout << "The path is: ";
            solutionGraph.Path.clear();
            solutionGraph.visit.resize(nVertex, false);
            solutionGraph.dfs(i);
            cout << '\n';
            counter ++;
            int cnt = 0;
            vector<int> temp_s, temp_t;
            for(const int &u: solutionGraph.Path){
                cout << u << ' ';
                if(u < nVertex / 2) {
                    temp_s.push_back(u);
                }else{
                    temp_t.push_back(u);
                }
                visited[u] = true;
            }
            assert(temp_s.size() % (K / 2) == 0);
            assert(temp_t.size() % (K / 2) == 0);
            sequenceSplit(temp_s, temp_t, kPath, K / 2, nVertex / 2);
            cout << '\n';
        }
    }
};

