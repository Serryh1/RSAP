#ifndef Input_h
#define Input_h

#include <utility>
#include <vector>
#include <string>
using std::vector;
using std::pair;
using std::string;



void ALG1_InputfileProcess(const string &str, int &nVertex, int &nEdges, int &nTerminals, 
    int *&terminals, vector<vector<double> > & adjMatrix, vector<vector<int >> & sites, vector<vector<int > > & E, int &VehicleNum, int &RequestNum, 
    vector<pair<double,double>> &Rs, vector<pair<double,double>> &Rt, vector<pair<double,double>> &VehilcePosition);

void ALG2_InputfileProcess(const string &str, int &nVertex, int &nEdges, int &nTerminals, 
    int *&terminals, vector<vector<double> > & adjMatrix, vector<vector<int >> & sites, vector<vector<int > > & E, int &VehicleNum, int &RequestNum, 
    vector<pair<double,double>> &Rs, vector<pair<double,double>> &Rt, vector<pair<double,double>> &VehilcePosition);

#endif