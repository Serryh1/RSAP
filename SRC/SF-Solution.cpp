#include "SF-Solution.h"
#include <unordered_map>
#include <set>
using namespace std;

steinerForest::steinerForest(int V)
{
    this->V = V;
    adj = new vector<node>[V];
}

void steinerForest::addEdge(int v, int w, double wei)
{
    node vNode(v,wei);
    node WNode(w, wei);
    adj[v].push_back(WNode); // Add w to v’s list.
    adj[w].push_back(vNode); // Add v to w’s list.
}

// delete this edge
void steinerForest::deleteEdge(int u, int v){
    int x = -1;
    for(int i = 0; i < adj[u].size(); ++ i){
        if(adj[u][i].index == v){
            x = i;
            break;
        }
    }
    if(x == -1) return ;
    adj[u].erase(adj[u].begin() + x);
    int y = -1;
    for(int i = 0; i < adj[v].size(); ++ i){
        if(adj[v][i].index == u){
            y = i;
            break;
        }
    }
    if(y == -1) return ;
    adj[v].erase(adj[v].begin() + y);
}

bool steinerForest::BFS(int s,int t)
{
    int u;
    bool *visited = new bool[V];
    for(int i = 0; i < V; i++)
    visited[i] = false;
    
    queue<int> myQueue;
    visited[s] = true;
    myQueue.push(s);
    
    while(!myQueue.empty())
    {
        u = myQueue.front();
        myQueue.pop();
        
        for (int i=0; i< adj[u].size(); i++)
        {
            if ((!visited[adj[u][i].index]) && (adj[u][i].needed))
            {
                visited[adj[u][i].index] = true;
                myQueue.push(adj[u][i].index);
                if (adj[u][i].index==t) return true;
            }
        }
    }
    return false;
}




void steinerForest::BFSCleanSolution(int s,int t, double &sum)
{
    // in this algorithm, it finds the path between each pair in the initial solution and zero out the path for the next pair
    sum=0;
    int u;
    bool *visited = new bool[V];
    int *pred = new int[V];
    int *We = new int[V];
    
    for(int i = 0; i < V; i++){
        visited[i] = false;
        pred[i]=-1;
        We[i]=0;
    }
    
    queue<int> myQueue;
    visited[s] = true;
    myQueue.push(s);
    
    while(!myQueue.empty())
    {
        u = myQueue.front();
        myQueue.pop();
        for (int i=0; i< adj[u].size(); i++)
        {
            if ((!visited[adj[u][i].index]))
            {
                visited[adj[u][i].index] = true;
                pred[adj[u][i].index]=u;
                We[adj[u][i].index]=adj[u][i].weight;
                myQueue.push(adj[u][i].index);
                if (adj[u][i].index==t) {
                    queue<int> empty;
                    swap(myQueue,empty);
                    break;
                }
            }
        }
    }
    vector<vector<int>> pathVec;
    
    int currentV=t;
    while (pred[currentV]!=-1) {
        sum+=We[currentV];
        for (int i=0; i< adj[currentV].size(); i++)
        {
            if (adj[currentV][i].index==pred[currentV]) {
                vector<int> tPath;
                tPath.push_back(currentV);
                tPath.push_back(adj[currentV][i].index);
                pathVec.push_back(tPath);
                adj[currentV][i].weight=0;
            }
        }
        for (int i=0; i< adj[pred[currentV]].size(); i++)
        {
            if (adj[pred[currentV]][i].index==currentV) {

                adj[pred[currentV]][i].weight=0;
            }
        }
        currentV=pred[currentV];
    }
}

bool steinerForest::isvalidcut(const int &u, const  int &v, const int &RequestNum, const int &k){
    queue<int> q;
    unordered_map<int,int> visited;
    q.push(u);
    bool flag = false;
    visited[v] = 1;
    set<int> se;
    int cnt = 0;
    while(!q.empty()){
        int t = q.front();
        q.pop();
        if(visited[t]) continue;
        visited[t] = true;
        cnt ++;
        se.insert(t);
        for(const auto  &NODE: adj[t]){
            int to = NODE.index;
            if(visited[to]) continue;
            q.push(to);
        }
    }
    if(cnt % k > 0) return false;

    // check each pair terminal is in the same component
    for(const auto &val: se){
        int val_t = Terminals[val];
        if(!se.count(val_t)) return false;
    }
    return true;
}

void steinerForest::dfs(int u){
    visit[u] = true;
    // cout << u << ' ';
    Path.push_back(u);
    for (int i = 0; i < adj[u].size(); i++){
        int v = adj[u][i].index;
        if(!visit[v]){
            dfs(v);
        }
    }
}