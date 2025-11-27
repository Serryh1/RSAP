#ifndef HRA_h
#define HRA_h
#include <vector>
#include <utility>

using namespace std;
class HRA
{
    private:
    public:

    vector<vector<double>> Adj;
    double GetEdgeVal(const pair<double,double> &xa, const pair<double,double> &xb);
    double PrimCompute(const vector<int> &vec, const vector<pair<double,double>> &pos, vector<pair<int, int>>& mst_edges);
    vector<pair<int,int>> CompleteGraphConstruct_k1(const vector<pair<double, double>> &Rs, const vector<pair<double, double>> &Rt,
        const vector<vector<int>> &group_s);
    vector<pair<int,int>> CompleteGraphConstruct_k2(const vector<pair<double, double>> &Rs, const vector<pair<double, double>> &Rt,
       const vector<vector<int>> &group_s, const int &m);

    vector<pair<int,int>> BipartiteGraphConstruct(const vector<pair<double, double>> &Rs, const vector<pair<double, double>> &Rt,
        const vector<vector<int>> &group_s, const vector<vector<int>> &discardedGroup);
    
    vector<vector<int>> GroupSplitResult(const int &k, const vector<pair<double, double>> &Rs, const vector<pair<double, double>> &Rt, const int &RequestNum, 
        const int &m);

    vector<int> VechileMatchResult(const vector<pair<double,double>> &VehilcePosition, const vector<vector<int>> &groups, const vector<pair<double,double>> &Rs);

    void HRASolve(const int &k, const vector<pair<double,double>> &VehilcePosition, const vector<pair<double, double>> &Rs, const vector<pair<double, double>> &Rt, const int &RequestNum, 
        const int &m);

};

#endif