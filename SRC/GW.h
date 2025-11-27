#ifndef GW_hpp
#define GW_hpp

#include <queue>

using namespace std;
void GW_ALG1(int &nVertex, int &nEdges, int &nTerminals, int *&terminals, vector<vector<double > > & adjMatrix, 
    vector<vector<int > > & pairedTerminals, vector<vector<int > > & E, const int &k, vector<vector<int>> &kPath);

void GW_ALG2(int &nVertex, int &nEdges, int &nTerminals, int *&terminals, vector<vector<double > > & adjMatrix, 
    vector<vector<int > > & pairedTerminals, vector<vector<int > > & E , const int &k, vector<vector<int>> &kPath);

#endif /* GW_hpp */
