#ifndef BipartiteGraph_h
#define BipartiteGraph_h

#include <vector>
#include <stack>
using namespace std;

class KM{
    public:
    vector<double> la, lb;
    vector<bool> va, vb;
    vector<int> match;
    vector<vector<double >> w;
    vector<double> slack;
    int n;
    double delta;
    void init(int n);
    bool dfs(int u);
    double solve();
};



void Alg1_bigraphCreate(const vector<pair<double,double>> &VehilcePosition, const vector<vector<int>> &kPath, 
    const vector<pair<double,double>> &Rs, const vector<pair<double,double>> &Rt, vector<vector<double>> &Adj);

void Alg2_bigraphCreate(const vector<pair<double,double>> &VehilcePosition, const vector<vector<int>> &kPath, 
    const vector<pair<double,double>> &Rs, const vector<pair<double,double>> &Rt, vector<vector<double>> &Adj, const int &gap);

void VehicleminCostMatching(const vector<vector<double>> &Adj, vector<int> &Result);

void GroupminCostMatching(const vector<vector<double>> &Adj, vector<int> &Result);

void Alg1_computeResult(const vector<pair<double,double>> &VehilcePosition, const vector<vector<int>> &kPath, const vector<pair<double,double>> &Rs, 
    const vector<pair<double,double>> &Rt, const vector<int> &Result);

void Alg2_computeResult(const vector<pair<double,double>> &VehilcePosition, const vector<vector<int>> &kPath, const vector<pair<double,double>> &Rs, 
    const vector<pair<double,double>> &Rt, const vector<int> &Result, const int &gap);
#endif