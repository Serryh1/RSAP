#ifndef SF_Solution_h
#define SF_Solution_h
#include <iostream>
#include <queue>
#include <unordered_map>

using namespace std;

class steinerForest
{
    
    public:
    steinerForest(int V);
    void addEdge(int v, int w,double wei);
    void deleteEdge(int u, int v);
    bool BFS(int s,int t);
    void BFSCleanSolution(int s,int t, double &sum);
    class node
    {
        public:
        node(int i, double w)
        {
            index=i;
            needed=true;
            weight=w;
        }
        int index;
        bool needed;
        double weight;
        
    };
    int V;    // No. of vertices, initialy we add all vertices
    vector<node> *adj; // Pointer to an array containing adjacency lists
    vector<int> Path;
    void dfs(int s);
    vector<bool> visit;
    unordered_map<int,int> Terminals;   
    bool isvalidcut(const int &u,const  int &v, const int &RequestNum, const int &k); // edge(u,v) can be cut
};

#endif /* SF_Solution_h */
