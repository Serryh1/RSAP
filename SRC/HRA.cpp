#include "HRA.h"
#include "BipartiteGraph.h"

#include <vector>
#include <cmath>
#include <utility>
#include <cassert>
#include <cfloat>
#include <algorithm>
#include <iostream>
#include <stack>
#include <cstring>
#include <queue>
#include <iomanip>
#include <unordered_map>
#include <set>
using namespace std;

#define pb push_back
#define mp make_pair
#define pii pair<int,int>
// #define M 4000
#define N 18010
#define se second
#define Base 100000000000
#define fi first
using T= long long;
const T INF=0x3f3f3f3f3f3f3f3f;

struct blossom_tree{

    struct edge{
        int u,v;
        T w;
        edge(){}
        edge(int u,int v,T w):u(u),v(v),w(w){}
    };
    // Graph
    int n,n_x; // [1, n]: point; [n+1, n_x]: flower
    // edge g[N][N]; // adjacent matrix
    vector<vector<edge>> g;  // using the vector to avoid MLE
    // flower
    vector<int>flower[N]; // nodes in flower i (outer flower)
    int root[N]; // flower root, root<=n root=i: normal nodes
    // int flower_from[N][N]; // flower_from[b][x]: outermost flower in b that contains x
    vector<vector<int>> flower_from;
    // slack
    T label[N*2]; // node label, [1, n] point label, [n+1, n_x] flower label
    int col[N]; // color saved at flower root
    int slv[N]; // slack node of NON LEAF NODES, slv[y]=x z(x,y) min_x
    // match
    int mat[N]; // match, mat[x]=y (x,y)\in E
    int fa[N]; // fa in cross tree
    int vis[N]; // if in path
    
    queue<int>Q; // bfs queue
    
    // calculate slv
    inline T calc_slv(edge e){return label[e.u]+label[e.v]-e.w;}
    inline void update_slv(int u,int v){if(!slv[v]||calc_slv(g[u][v])<calc_slv(g[slv[v]][v]))slv[v]=u;}
    inline void recalc_slv(int u){
        slv[u]=0;
        for(int i=1;i<=n;i++)if(g[i][u].w>0&&root[i]!=u&&col[root[i]]==1)update_slv(i,u);
    }
    
    // only push nodes, not flowers
    void q_push(int x){
        if(x<=n)Q.push(x);
        else for(auto p:flower[x])q_push(p);
    }
    
    // set root of all nodes in x to r
    void set_root(int x,int r){
        root[x]=r;
        if(x>n)for(auto p:flower[x])set_root(p,r);
    }
    
    // return a (+-)^k path in flower b from root[b] to x
    int get_even_path_in_flower(int b,int x){
        int pr=find(flower[b].begin(),flower[b].end(),x)-flower[b].begin();
        // b is flower, x in b
        if(pr%2==0)return pr;
        reverse(flower[b].begin()+1,flower[b].end());
        return flower[b].size()-pr;
    }
    
    // set (u,v) match, can be flower
    void set_match(int u,int v){
        mat[u]=g[u][v].v;
        if(u>n){
            edge e=g[u][v];
            int xr=flower_from[u][e.u];
            int pr=get_even_path_in_flower(u,xr);
            for(int i=0;i<pr;i++)set_match(flower[u][i],flower[u][i^1]);
            set_match(xr,v);
            rotate(flower[u].begin(),flower[u].begin()+pr,flower[u].end()); // change receptacle
        }
    }
    
    // link 2 S points
    void side_augment(int u,int v){
        int nv=root[mat[u]],nu=root[fa[nv]];
        while(1){
            set_match(u,v);
            u=nu,v=nv;
            if(!nv)break;
            set_match(nv,nu);
            nv=root[mat[u]],nu=root[fa[nv]];
        }
    }
    void linkSS(int u,int v){
        side_augment(u,v); 
        side_augment(v,u);
    }
    
    int get_lca(int u,int v){
        static int t=0;
        ++t; // to avoid clearing vis
        while(u||v){
            if(vis[u]==t)return u;
            vis[u]=t;
            u=root[mat[u]];
            if(u)u=root[fa[u]];
            if(!u)swap(u,v);
        }
        return 0;
    }
    
    void add_blossom(int u,int v,int r){
        int i,b=n+1;
        while(b<=n_x&&root[b])b++;
        if(b>n_x)++n_x;
        // clear
        col[b]=1;label[b]=0;mat[b]=mat[r];flower[b].clear();
        for(i=1;i<=n_x;i++)g[i][b].w=g[b][i].w=0;
        for(i=1;i<=n;i++)flower_from[b][i]=0;
        // construct flower
        while(u!=r){
            flower[b].pb(u);u=root[mat[u]];q_push(u);
            flower[b].pb(u);u=root[fa[u]];
        }
        flower[b].pb(r);
        reverse(flower[b].begin(),flower[b].end());
        while(v!=r){
            flower[b].pb(v);v=root[mat[v]];q_push(v);
            flower[b].pb(v);v=root[fa[v]];
        }
        // set as outermost flower
        set_root(b,b);
        // calculate slack
        for(auto p:flower[b]){
            for(i=1;i<=n_x;i++){
                // set to min slave
                if(!g[b][i].w||calc_slv(g[p][i])<calc_slv(g[b][i])){
                    g[b][i]=g[p][i];
                    g[i][b]=g[i][p];
                }
            }
            for(i=1;i<=n;i++)if(flower_from[p][i])flower_from[b][i]=p;
        }
        recalc_slv(b);
    }
    
    // only expand outermost blossom b, b is T(white) blossom
    void expand_blossom(int b){
        int i,x;
        for(auto p:flower[b])set_root(p,p);
        x=flower_from[b][g[b][fa[b]].u];
        // [0,pr]: (+-)^k, insert into tree, add black to queue
        int pr=get_even_path_in_flower(b,x);
        col[x]=2;fa[x]=fa[b];
        for(i=0;i<pr;i+=2){
            // from bottom to upper layer in tree
            int white=flower[b][i];
            int black=flower[b][i+1];
            col[black]=1;col[white]=2;
            fa[white]=g[black][white].u;
            slv[black]=slv[white]=0;
            q_push(black);
        }
        // others: color=0
        for(i=pr+1;i<flower[b].size();i++){
            col[flower[b][i]]=0;
            recalc_slv(flower[b][i]);
        }
        // delete b
        root[b]=0;
        flower[b].clear();
    }
    
    // found_edge
    int augment_path(edge e){
        int u=root[e.u],v=root[e.v];
        if(!col[v]){
            fa[v]=e.u;
            col[v]=2;
            int nu=root[mat[v]];
            slv[nu]=slv[v]=0;
            col[nu]=1;
            q_push(nu);
        }else if(col[v]==1){
            int r=get_lca(u,v);
            if(r)add_blossom(u,v,r);
            else return linkSS(u,v),1;
        }
        return 0;
    }
    
    int augment(){
        int i;
        memset(col,0,sizeof(int)*(n_x+1));
        memset(slv,0,sizeof(int)*(n_x+1));
        memset(fa,0,sizeof(int)*(n_x+1));
        Q=queue<int>();
        for(i=1;i<=n_x;i++)
            if(root[i]==i&&!mat[i]){
                // add all unmatched points
                col[i]=1;
                q_push(i);
            }
        if(Q.empty())return 0;
        while(1){
            while(!Q.empty()){
                int p=Q.front();Q.pop();
                for(i=1;i<=n;i++){
                    if(g[p][i].w==0||root[i]==root[p])continue;
                    // not in same flower
                    T d=calc_slv(g[p][i]);
                    if(!d){if(augment_path(g[p][i]))return 1;}
                    else if(col[root[i]]!=2)update_slv(p,root[i]);
                }
            }
            T delta=INF;
            // calc delta
            for(i=1;i<=n;i++)if(col[root[i]]==1)delta=min(delta,label[i]);
            for(i=n+1;i<=n_x;i++)if(root[i]==i&&col[i]==2)delta=min(delta,label[i]/2);
            for(i=1;i<=n_x;i++){
                if(root[i]!=i||!slv[i])continue;
                if(!col[i])delta=min(delta,calc_slv(g[slv[i]][i]));
                else if(col[i]==1)delta=min(delta,calc_slv(g[slv[i]][i])/2);
            }
            // update label
            for(i=1;i<=n;i++){
                if(col[root[i]]==1)label[i]-=delta;
                else if(col[root[i]]==2)label[i]+=delta;
            }
            for(i=n+1;i<=n_x;i++){
                if(root[i]!=i)continue;
                if(col[i]==1)label[i]+=2*delta;
                else if(col[i]==2)label[i]-=2*delta;
            }
            for(i=1;i<=n;i++)if(label[i]<=0)return 0;
            for(i=1;i<=n_x;i++){
                if(root[i]!=i||!slv[i]||root[slv[i]]==i)continue;
                if(calc_slv(g[slv[i]][i])==0&&augment_path(g[slv[i]][i]))return 1;
            }
            // expand
            for(i=n+1;i<=n_x;i++)
                if(root[i]==i&&col[i]==2&&label[i]==0)
                    expand_blossom(i);
        }
        return 0;
    }
    
    void init(int _n,vector<pair<T,pii>>edges){
        int i,j;
        n=n_x=_n;
        memset(mat,0,sizeof(mat));
        for(i=0;i<=n;i++){
            root[i]=i;
            flower[i].clear();
            for(j=0;j<=n;j++){
                flower_from[i][j]=(i==j)?i:0;
                g[i][j]=edge(i,j,0);
            }
        }
        T w_max=-1e9;
        for(auto pr:edges){
            int u=pr.se.fi,v=pr.se.se;
            T w=pr.fi;
            g[u][v]=edge(u,v,w*2);
            g[v][u]=edge(v,u,w*2);
            w_max=max(w_max,w);
        }
        for(i=1;i<=n;i++)label[i]=w_max;
    }
    
    pair<int,T>calc(){
        int i,cnt=0;T s=0;
        while(augment())++cnt;
        for(i=1;i<=n;i++)if(mat[i]>i)s+=g[i][mat[i]].w/2;
        return mp(cnt,s);
    }
    void set_up(int n){
        g.resize(n, vector<edge>(n));
        flower_from.resize(n, vector<int>(n));
    }
};

double HRA::GetEdgeVal(const pair<double,double> &xa, const pair<double,double> &xb){
    double dist = sqrt((xa.first - xb.first) * (xa.first - xb.first) + (xa.second - xb.second) * (xa.second - xb.second));
    return dist;
}

double HRA::PrimCompute(const vector<int> &vec, const vector<pair<double,double>> &pos, vector<pair<int, int>>& mst_edges){

    int n = vec.size();
    vector<vector<double>> cost(n, vector<double>(n));

    for(int i = 0; i < n; ++ i){
        for(int j = 0; j < n; ++ j){
            int val_i = vec[i];
            int val_j = vec[j];
            double edgecost = GetEdgeVal(pos[val_i], pos[val_j]);
            cost[i][j] = edgecost;
        }
    }

    vector<int> parent(n, -1);
    vector<double> key(n, DBL_MAX);  
    vector<bool> inMST(n, false);

    key[0] = 0;  

    for (int count = 0; count < n - 1; ++count) {
        int u = -1;
        double min_key = DBL_MAX;
        for (int v = 0; v < n; ++v) {
            if (!inMST[v] && key[v] < min_key) {
                min_key = key[v];
                u = v;
            }
        }

        if (u == -1) {
            cerr << "Error: Graph is disconnected, no MST exists." << endl;
            mst_edges.clear();  
            return -1.0;
        }

        inMST[u] = true;

        for (int v = 0; v < n; ++v) {
            if (!inMST[v] && cost[u][v] > 0 && cost[u][v] < key[v]) {
                key[v] = cost[u][v]; 
                parent[v] = u;         
            }
        }
    }

    mst_edges.clear();
    double total_weight = 0.0;
    for (int i = 1; i < n; ++i) {  
        if (parent[i] != -1) {
            mst_edges.emplace_back(parent[i], i);  // addedge (parent[i], i)
            total_weight += cost[parent[i]][i];   // 
        }
    }
    return total_weight;
}
blossom_tree blossom;

vector<pair<int,int>> HRA::CompleteGraphConstruct_k1(const vector<pair<double, double>> &Rs, const vector<pair<double, double>> &Rt,const vector<vector<int>> &group_s){
    int sz = group_s.size();

    int siz_group = Adj[0].size();
    
    vector<vector<double>> ADJ(sz, vector<double>(sz));

    for(int i = 0; i < sz; ++ i){
        auto vecs_i = group_s[i];
        for(int j = i + 1; j < sz; ++ j){
            auto vecs_j = group_s[j];
            auto temp = vecs_i;
            temp.insert(temp.end(), vecs_j.begin(), vecs_j.end());
            vector<pair<int, int>> mst_edges;
            double MST_s = PrimCompute(temp, Rs, mst_edges);
            double MST_t = PrimCompute(temp, Rt, mst_edges);
            double MST_st = MST_s + MST_t;

            double MST_i_st = PrimCompute(vecs_i, Rs, mst_edges) + PrimCompute(vecs_i, Rt, mst_edges);
            double MST_j_st = PrimCompute(vecs_j, Rs, mst_edges) + PrimCompute(vecs_j, Rt, mst_edges);

            MST_st = MST_st - MST_i_st - MST_j_st;

            int x = ceil(MST_st);
            ADJ[i][j] = Base - x;
            ADJ[j][i] = Base - x;
        }
    }

    vector<pair<int,int>> result;


    vector<pair<T,pii>>edges;

    // 
    for(int i = 0; i < sz ; i ++){
        for(int j = i+1; j < sz ; j ++){
            edges.push_back(make_pair((T)ADJ[i][j], make_pair(i + 1, j + 1)));
        }
    }

    blossom.init(sz, edges);
    double ans = blossom.calc().second;
    for(int i = 1; i <= sz; ++ i){
        if(i < blossom.mat[i]){
            result.push_back(make_pair(i - 1, blossom.mat[i] - 1));
        }
        cout << i << ' ' << blossom.mat[i] << endl;
    }

    // printf("beautiful match's cost %.6lf\n", Base * ((sz) / 2) - ans);
    cout << "beautiful match's cost  is: " <<  fixed << setprecision(6) << Base * ((sz) / 2) - ans << endl;
    return result;

}

vector<pair<int,int>> HRA::CompleteGraphConstruct_k2(const vector<pair<double, double>> &Rs, const vector<pair<double, double>> &Rt,
        const vector<vector<int>> &group_s, const int &m){
    int sz = group_s.size();
    int siz_group = Adj[0].size();

    vector<vector<double>> ADJ(sz+m, vector<double>(sz + m));

    for(int i = 0; i < sz; ++ i){
        auto vecs_i = group_s[i];
        for(int j = i + 1; j < sz; ++ j){
            auto vecs_j = group_s[j];
            auto temp = vecs_i;
            temp.insert(temp.end(), vecs_j.begin(), vecs_j.end());
            vector<pair<int, int>> mst_edges;
            double MST_s = PrimCompute(temp, Rs, mst_edges);
            double MST_t = PrimCompute(temp, Rt, mst_edges);
            double MST_st = MST_s + MST_t;

            double MST_i_st = PrimCompute(vecs_i, Rs, mst_edges) + PrimCompute(vecs_i, Rt, mst_edges);
            double MST_j_st = PrimCompute(vecs_j, Rs, mst_edges) + PrimCompute(vecs_j, Rt, mst_edges);

            MST_st = MST_st - MST_i_st - MST_j_st;

            int x = ceil(MST_st);
            ADJ[i][j] = Base - x;
            ADJ[j][i] = Base - x;
        }
    }

    // we add another m vertices(sz+1, ..., sz + m)
    // addedge(sz + i, sz + j, INF), addedge(i, sz + j, 0); 
    // that is we must choose m pairwise(i, sz + j)  and not choose(sz+i, sz+j)
    // we use maximum value matching, so the value should be inf and 0

    for(int i = sz; i < sz + m; ++ i){
        for(int j = sz; j < sz + m; ++ j){
            ADJ[i][j] = 0;
            ADJ[j][i] = 0;
        }
        for(int j = 0; j < sz; ++ j){
            ADJ[i][j] = Base;
            ADJ[j][i] = Base;
        }
    }

    vector<pair<int,int>> result;
    // blossom_tree blossom;

    vector<pair<T,pii>>edges;
    // 
    for(int i = 0; i < sz + m; i ++){
        for(int j = i+1; j < sz + m; j ++){
            edges.push_back(make_pair((T)ADJ[i][j], make_pair(i + 1, j + 1)));
        }
    }

    blossom.init(sz + m, edges);
    double ans = blossom.calc().second;
    int cnt = 0;
    for(int i = 1; i <= sz + m; ++ i){
        if(i < blossom.mat[i]){
            result.push_back(make_pair(i - 1, blossom.mat[i] - 1));
        }
        cout << i << ' ' << blossom.mat[i] << endl;
    }
    // printf("beautiful match's cost %.6lf\n", Base * ((sz + m) / 2) - ans);
    cout << "beautiful match's cost  is: " <<  fixed << setprecision(6) << Base * ((sz + m) / 2) - ans << endl;
    return result;
}

vector<pair<int,int>> HRA::BipartiteGraphConstruct(const vector<pair<double, double>> &Rs, const vector<pair<double, double>> &Rt,
        const vector<vector<int>> &group_s, const vector<vector<int>> &discardedGroup){
    
    int sz = group_s.size();

    vector<vector<double>> ADJ(sz, vector<double>(sz));

    for(int i = 0; i < sz; ++ i){
        vector<int> vecs_i = group_s[i];
        for(int j = i; j < sz; ++ j){
            vector<int> vecs_j = discardedGroup[j];
            auto temp = vecs_i;
            temp.insert(temp.end(), vecs_j.begin(), vecs_j.end());
            vector<pair<int, int>> mst_edges;
            double MST_s = PrimCompute(temp, Rs, mst_edges);
            double MST_t = PrimCompute(temp, Rt, mst_edges);
            double MST_st = MST_s + MST_t;

            double MST_i_st = PrimCompute(vecs_i, Rs, mst_edges) + PrimCompute(vecs_i, Rt, mst_edges);
            double MST_j_st = PrimCompute(vecs_j, Rs, mst_edges) + PrimCompute(vecs_j, Rt, mst_edges);

            MST_st = MST_st - MST_i_st - MST_j_st;

            ADJ[i][j] = MST_st;
            ADJ[j][i] = MST_st;
        }
    }

    vector<pair<int,int>> result;
    vector<int> temp;
    GroupminCostMatching(ADJ, temp);

    for(int i = 0; i < temp.size(); ++ i){
        result.push_back(make_pair(temp[i], i));
    }
    return result;
} 

vector<vector<int>> HRA::GroupSplitResult(const int &k, const vector<pair<double, double>> &Rs, const vector<pair<double, double>> &Rt, const int &mk, const int &m){
    

    cout << "GroupSpliting....." << endl;
    std::vector<int> binarycode;
    int num = k;
    int cnt1 = 0;
    while(num){
        if(num % 2 == 1) cnt1 ++ ;
        binarycode.push_back(num % 2);
        num /= 2;
    }
    // reverse(binarycode.begin(), binarycode.end());
    // that is k = 2^l
    vector<vector<int>> group_s;
    
    int l = binarycode.size() ;
    if(cnt1 == 1){

        group_s.resize(mk);
        for(int i = 0; i < mk; ++ i){ 
            vector<int> temp;
            temp.push_back(i);
            group_s[i]=temp;
        }

        // l-1 round combination
        for(int i = 0; i < l - 1; ++ i){
            vector<pair<int,int>> matchResult;
            matchResult = CompleteGraphConstruct_k1(Rs, Rt, group_s);

            int now_groupsize = matchResult.size();

            vector<vector<int>> temp_group_s;
            
            for(const auto &[x, y]: matchResult){
                vector<int> temp;
                for(const auto &v: group_s[x]){
                    temp.push_back(v);
                }
                for(const auto &v: group_s[y]){
                    temp.push_back(v);
                }
                temp_group_s.push_back(temp);
            }

            // combination complete
            // cover origin group_s
            group_s = temp_group_s;
        }

        // this time, each group in group_s, its size is k
        return group_s;
    }else{
        
        group_s.resize(mk);
        for(int i = 0; i < mk; ++ i){ 
            vector<int> temp;
            temp.push_back(i);
            group_s[i]=temp;
        }

        stack<vector<vector<int>>> discardedGroups;
        // 1101 -----> 1011  reverse
        // split process
        for(int i = 0; i < l - 1; ++ i){
            // that is this time we will discard m groups
            if(binarycode[i]){
                vector<pair<int,int>> matchResult;
                int discardIdx = group_s.size();
                matchResult = CompleteGraphConstruct_k2(Rs, Rt, group_s, m);

                int match_groupsize = matchResult.size();
                vector<vector<int>> temp_group_s;
                vector<vector<int>> discard_group_s;

                for(const auto &[x, y]: matchResult){
                    vector<int> temp;
                    if(x < discardIdx && y < discardIdx){    
                        for(const auto &v: group_s[x]){
                            temp.push_back(v);
                        }
                        for(const auto &v: group_s[y]){
                            temp.push_back(v);
                        }
                        temp_group_s.push_back(temp);
                    }else{
                        int v = min(x, y);
                        for(const auto &v: group_s[v]){
                            temp.push_back(v);
                        }
                        discard_group_s.push_back(temp);
                    }
                }

                discardedGroups.push(discard_group_s);
                // update now group_s
                group_s = temp_group_s;
            }else{
                // 
                vector<pair<int,int>> matchResult;
                matchResult = CompleteGraphConstruct_k2(Rs, Rt, group_s, 0);

                int match_groupsize = matchResult.size();

                vector<vector<int>> temp_group_s;
                
                for(const auto &[x, y]: matchResult){
                    vector<int> temp;
                    for(const auto &v: group_s[x]){
                        temp.push_back(v);
                    }
                    for(const auto &v: group_s[y]){
                        temp.push_back(v);
                    }
                    temp_group_s.push_back(temp);
                }
                group_s = temp_group_s;
            }   
        }
        // combination process
         
        while(!discardedGroups.empty()){
            auto vec = discardedGroups.top();
            discardedGroups.pop();

            auto matchResult = BipartiteGraphConstruct(Rs, Rt, group_s, vec);
            int match_groupsize = matchResult.size();

            vector<vector<int>> temp_group_s;
            for(const auto &[x, y]: matchResult){
                vector<int> temp;
                for(const auto &v: group_s[x]){
                    temp.push_back(v);
                }
                for(const auto &v: vec[y]){
                    temp.push_back(v);
                }
                temp_group_s.push_back(temp);
            }
            group_s = temp_group_s;
        }
        return group_s;
    }
    
}


vector<int> HRA::VechileMatchResult(const vector<pair<double,double>> &VehilcePosition, const vector<vector<int>> &groups, 
    const vector<pair<double,double>> &Rs){

    int VehicleSize = VehilcePosition.size();
    int GroupSize = groups.size();

    vector<vector<double>> cost(VehicleSize, vector<double> (VehicleSize)); 
    
    for(int i = 0; i < VehicleSize; ++ i){
        auto vpos = VehilcePosition[i];
        for(int j = 0; j < GroupSize; ++ j){
            auto group_j = groups[j];
            double MinDist = DBL_MAX;
            for(int u: group_j){
                auto pos_u = Rs[u];
                double dist = GetEdgeVal(vpos, pos_u);
                if(dist - MinDist <= 1e-4){
                    MinDist = dist;
                }
            }
            // Vehicle -- group
            cost[i][j] = MinDist;
        }
    }

    vector<int> MatchResult;

    VehicleminCostMatching(cost, MatchResult);
    return MatchResult;
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

void HRA::HRASolve(const int &k, const vector<pair<double,double>> &VehilcePosition, const vector<pair<double, double>> &Rs, 
    const vector<pair<double, double>> &Rt, const int &RequestNum, const int &m){
    
    vector<vector<int>> groups;
    blossom.set_up(2*(RequestNum + m + 1));
    

    Adj.resize(RequestNum + m, vector<double>(RequestNum + m));
    // m is the cehicle number
    groups = GroupSplitResult(k, Rs, Rt, RequestNum, m);

    vector<int> vehiclematch = VechileMatchResult(VehilcePosition, groups, Rs);

    unordered_map<int,int> mp;
    for(int i = 0; i < vehiclematch.size(); ++ i){
        cout << "Vehicle number is :" << vehiclematch[i] << " --- group id is: "<< i << '\n';
        mp[vehiclematch[i]] = i;
    }

    for(int i = 0; i < m; ++ i){
        cout << "group id: " << i << " The requests are: ";
        for(int j : groups[i]){
            cout << j << ' ';
        }
        cout << '\n';
    }

    // For each vehicle -- group
    // we design the path
    double tot_cost = 0;
    for(int VechicleId = 0; VechicleId < m; ++ VechicleId){
        int Start_S = -1;
        int End_S = -1;
        int Connect_T = -1;
        auto vpos = VehilcePosition[VechicleId];
        double MinDist = INF;

        auto vec = groups[mp[VechicleId]];

        for(const auto &v: vec){
            double dist = GetEdgeVal(vpos, Rs[v]);
            if(dist - MinDist <= 1e-4){
                MinDist = dist;
                Start_S = v;
            }
        }
        double tot = 0;
        // u -- Start_S
        tot = 0;
        tot += MinDist;

        int sz = vec.size();
        MinDist = INF;
        for(int i = 0; i < sz; ++ i){
            for(int j = 0; j < sz; ++ j){
                double dist = GetEdgeVal(Rs[vec[i]], Rt[vec[j]]);
                if(dist - MinDist <= 1e-4){
                    MinDist = dist;
                    End_S = vec[i];
                    Connect_T = vec[j];
                }
            }
        }

        cout << "Start_S is :" << Start_S << endl;
        cout << "End_S is :" << End_S << endl;
        cout << "Connect_T is :" << Connect_T << endl;
        
        vector<int> Path;
        // V 
        Path.push_back(VechicleId);

        vector<pair<int, int>> mst_s_edges;
        double Mst_s = PrimCompute(vec, Rs, mst_s_edges);

        // Graph g_s(RequestNum);
        std::unordered_map<int, std::vector<int>> g_s;
        for(const auto &[u, v]: mst_s_edges){
            g_s[vec[u]].push_back(vec[v]);
            g_s[vec[v]].push_back(vec[u]);
        }


        auto getpath = [&](int start, int end, int max_idx, std::unordered_map<int, std::vector<int>> g_s){
            vector<int> pre(max_idx, -1);
            vector<bool> visted(max_idx, 0);
            queue<int> q;
            q.push(start);
            while(!q.empty()){
                int t = q.front();
                q.pop();
                if(visted[t]) continue;
                visted[t] = true;
                if(t == end) return pre;
                for(auto &to: g_s[t]){
                    if(visted[to]) continue;
                    pre[to] = t;
                    q.push(to);
                }
            }
            return pre;
        };

        vector<int> pre;

        // find the shortest path from Start_S to End_S
        pre = getpath(Start_S, End_S, RequestNum, g_s);

        vector<int> path;
        int ed = End_S;
        path.push_back(ed);
        while(pre[ed] != -1){
            ed = pre[ed];
            path.push_back(ed);
        }
        reverse(path.begin(), path.end());
        
        auto func = [&](int x, int y){
            for(int i = 1; i < path.size(); ++ i){
                if(path[i-1] == x && path[i] == y) {
                    return false;
                }
                if(path[i-1] == y && path[i] == x) {
                    return false;
                }
            }
            return true;
        };
        // double the other edges escape the shortest path
        for(const auto &[u, v]: mst_s_edges){
            if(func(vec[u], vec[v])){
                // g_s.addEdge(vec[u], vec[v], 0);
                g_s[vec[u]].push_back(vec[v]);
                g_s[vec[v]].push_back(vec[u]);
            }
        }
        // Test
        // g_s.printGraph();
        vector<int> Spath;
        if (1) {
            cout << "\nEulerian Path Exists!" << endl;             
            vector<int> temp;
            // solve the path from Start_S 
            dfsEuler(Start_S, g_s, temp);
            cout << "Path: ";
            int last = -1;
            for (int node : temp) {
                cout << node << " ";
            }
            cout << endl;
            set<int> se;

            cout << "Remove the same nodes:" << endl;
            pair<double,double> last_pos;
            bool ff = false;
            for(const auto &v: temp){
                if(se.count(v)) continue;
                se.insert(v);
                Spath.push_back(v);
                int op, idx;
                pair<double,double> now;
                cout << "S" << v << "->";
                idx = v;
                now = Rs[idx];
                // check whether is the start point
                if(ff) {
                    tot += GetEdgeVal(last_pos, now);
                }
                last_pos = now;
                ff = true;
            }
            // end_s ---> Connect_T
            tot += GetEdgeVal(last_pos, Rt[Connect_T]);
            cout << endl;
            cout << "S Path's Cost is :" << tot << endl;
        }else {
            cout << "\nNo Eulerian Path Exists!" << endl;
        }
        
       std::unordered_map<int, std::vector<int>> g_t;
        vector<pair<int, int>> mst_t_edges;
        double Mst_t = PrimCompute(vec, Rt, mst_t_edges);

        for(const auto &[u, v]: mst_t_edges){
            // g_t.addEdge(vec[u], vec[v], 0);
            // g_t.addEdge(vec[v], vec[u], 0);
            g_t[vec[u]].push_back(vec[v]);
            g_t[vec[v]].push_back(vec[u]);

        }
        // Test 
        // g_t.printGraph();
        vector<int> Tpath;
        if (1) {
            cout << "Eulerian Path Exists!" << endl;             
            // auto temp = g_t.findEulerPath(Connect_T);
            vector<int> temp;
            dfsEuler(Connect_T, g_t, temp);
            cout << "Path: ";
            int last = -1;
            for (int node : temp) {
                cout << node << " ";
            }
            cout << endl;
            set<int> se;
            cout << "Remove the same nodes:" << endl;
            pair<double,double> last_pos;
            bool ff = false;
            for(const auto &v: temp){
                if(se.count(v)) continue;
                se.insert(v);
                Tpath.push_back(v);
                int op, idx;
                pair<double,double> now;
                cout << "T" << v << "->";
                idx = v;
                now = Rt[idx];
                if(ff) {
                    tot += GetEdgeVal(last_pos, now);
                }
                last_pos = now;
                ff = true;
            }
            cout << endl;
            cout << "T Path's Cost is :" << tot << endl;
        }else {
            cout << "\nNo Eulerian Path Exists!" << endl;
        }

        cout << "The final Path is: \n";
        cout << "V" << VechicleId << "->";
        for(const auto &v: Spath){
            cout << "S" << v << "->";
        }
        for(const auto &v: Tpath){
            cout << "T" << v << "->";
        }
        cout << "\n";

        cout << "The cost for Vehicle "<< VechicleId << ", The group id is: " << mp[VechicleId] << endl;
        int last = Spath.size();
        tot_cost += tot;
        
    }
    cout << "---------------------\n";
    // cout << "TOT Cost is: " << tot_cost << setprecision(5) << num << endl;
    // printf("TOT COST is: %.6lf\n", tot_cost);
    cout << "TOT COST is: " << fixed << setprecision(6) << tot_cost << endl;
    cout << "----------------------\n";
    
}


