#include "Input.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <utility>
#include <cmath>

using namespace std;



void ALG1_InputfileProcess(const string &str, int &nVertex, int &nEdges, int &nTerminals, int *&terminals, vector<vector<double> > & adjMatrix, 
    vector<vector<int >> & sites, vector<vector<int > > & E, int &VehicleNum, int &RequestNum, vector<pair<double,double>> &Rs, 
    vector<pair<double,double>> &Rt, vector<pair<double,double>> &VehilcePosition){

    sites.clear();
    ifstream myfile;
    myfile.open(str);
    if(!myfile)
    {
        cout << "Error: file could not be opened"  << endl;
        exit (EXIT_FAILURE);
    }
    string s;
    const int &max_line = 65536;

    myfile >> s; // Vehicle
    // cout << s << ' ';
    myfile >> VehicleNum;
    // cout << VehicleNum << '\n';

    for(int i = 0; i < VehicleNum; ++ i){
        myfile >> s; // v
        string x, y;
        myfile >> x >> y;
        // cout << x << ' ' << y << '\n';
        double index_x = stod(x);
        double index_y = stod(y);
        VehilcePosition.push_back(make_pair(index_x, index_y)); // each vehicle's position
    }
    
    myfile >> s; // Request
    // cout << s << '\n';
    myfile >> RequestNum; 

    nVertex = RequestNum;  // Graph H ---> nVertex
    nEdges = (nVertex - 1) * nVertex / 2; // H is a complete graph
    vector<double> tempAdjVec (nVertex);
    for (int i=0; i<nVertex; i++)
        adjMatrix.push_back(tempAdjVec);
    // Request = k*VechileNum;

    // cout << "s \n";
    for(int i = 0; i < RequestNum; ++ i){
        myfile >> s;  // s
        string x, y;
        myfile >> x >> y;
        // cout << x << ' ' << y << '\n';
        double index_x = stod(x);
        double index_y = stod(y);
        Rs.push_back(make_pair(index_x, index_y));
    }
    // cout << "t \n";
    for(int i = 0; i < RequestNum; ++ i){
        myfile >> s;  // t
        string x, y;
        myfile >> x >> y;
        // cout << x << ' ' << y << '\n';
        double index_x = stod(x);
        double index_y = stod(y);
        Rt.push_back(make_pair(index_x, index_y));
    }

    auto GetEdgeVal = [&](const int &i, const int &j){
        double Ws = (Rs[i].first - Rs[j].first)*(Rs[i].first - Rs[j].first) + (Rs[i].second - Rs[j].second)*(Rs[i].second - Rs[j].second);
        double WS = sqrt(Ws);
        double Wt = (Rt[i].first - Rt[j].first)*(Rt[i].first - Rt[j].first) + (Rt[i].second - Rt[j].second)*(Rt[i].second - Rt[j].second);
        double WT = sqrt(Wt) * 3 / 4;
        return WS + WT;
    };


    for(int i = 0; i < RequestNum; ++ i){
        for(int j = i + 1; j < RequestNum; ++ j){
            double dist = GetEdgeVal(i, j);
            adjMatrix[i][j] = dist;
            adjMatrix[j][i] = dist;
            vector<int> edge;
            edge.push_back(i);
            edge.push_back(j);
            E.push_back(edge); // edge set
        }
    }

    nTerminals = nVertex;
    terminals = new int[nVertex];

    for(int i = 0; i < nVertex; ++ i){
        terminals[i] = i;
        vector<int> terComponent;
	    terComponent.push_back(i);
	    terComponent.push_back(i);
	    sites.push_back(terComponent);
    }
}

// Alg2 input
void ALG2_InputfileProcess(const string &str, int &nVertex, int &nEdges, int &nTerminals, int *&terminals, vector<vector<double> > & adjMatrix, 
    vector<vector<int >> & sites, vector<vector<int > > & E, int &VehicleNum, int &RequestNum, vector<pair<double,double>> &Rs, 
    vector<pair<double,double>> &Rt, vector<pair<double,double>> &VehilcePosition){

    sites.clear();
    ifstream myfile;
    myfile.open(str);
    if(!myfile)
    {
        cout << "Error: file could not be opened"  << endl;
        exit (EXIT_FAILURE);
    }
    string s;
    const int &max_line = 65536;

    myfile >> s; // Vehicle
    myfile >> VehicleNum;

    for(int i = 0; i < VehicleNum; ++ i){
        myfile >> s; // v
        string x, y;
        myfile >> x >> y;
        double index_x = stod(x);
        double index_y = stod(y);
        VehilcePosition.push_back(make_pair(index_x, index_y)); // each vehicle's position
    }
    
    myfile >> s; // Request
    myfile >> RequestNum; 
    
    nVertex = RequestNum*2;  // Graph H ---> nVertex
    nEdges = (nVertex - 1) * nVertex / 2; // H is a complete graph
    vector<double> tempAdjVec (nVertex);
    for (int i=0; i<nVertex; i++)
        adjMatrix.push_back(tempAdjVec);
    // Request = k*VechileNum;

    // cout << "s \n";
    for(int i = 0; i < RequestNum; ++ i){
        myfile >> s;  // s
        string x, y;
        myfile >> x >> y;
        double index_x = stod(x);
        double index_y = stod(y);
        Rs.push_back(make_pair(index_x, index_y));
    }
    // cout << "t \n";
    for(int i = 0; i < RequestNum; ++ i){
        myfile >> s;  // t
        string x, y;
        myfile >> x >> y;
        double index_x = stod(x);
        double index_y = stod(y);
        Rt.push_back(make_pair(index_x, index_y));
    }

    auto GetEdgeVal = [&](const pair<double,double> &xa, const pair<double,double> &xb){
        double dist = sqrt((xa.first - xb.first) * (xa.first - xb.first) + (xa.second - xb.second) * (xa.second - xb.second));
        return dist;
    };


    for(int i = 0; i < RequestNum; ++ i){
        // addedge(si, sj)
        for(int j = i + 1; j < RequestNum; ++ j){
            double dis = GetEdgeVal(Rs[i], Rs[j]);
            adjMatrix[i][j] = dis;
            adjMatrix[j][i] = dis;
            vector<int> edge;
            edge.push_back(i);
            edge.push_back(j);
            E.push_back(edge); // edge set
        }
        // addedge(si, tj)
        // tj = sj + RequestNum
        for(int j = 0; j < RequestNum; ++ j){
            double dis = GetEdgeVal(Rs[i], Rt[j]);
            adjMatrix[i][j + RequestNum] = dis;
            adjMatrix[j + RequestNum][i] = dis;
            vector<int> edge;
            edge.push_back(i);
            edge.push_back(j + RequestNum);
            E.push_back(edge); // edge set
        }
        // addedge(ti, tj)
        for(int j = i + 1; j < RequestNum; ++ j){
            double dis = GetEdgeVal(Rt[i], Rt[j]);
            adjMatrix[i + RequestNum][j + RequestNum] = dis;
            adjMatrix[j + RequestNum][i + RequestNum] = dis;
            vector<int> edge;
            edge.push_back(i + RequestNum);
            edge.push_back(j + RequestNum);
            E.push_back(edge); // edge set
        }
    }

    nTerminals = RequestNum;
    terminals = new int[RequestNum];

    for(int i = 0; i < RequestNum; ++ i){
        terminals[i] = i;
        vector<int> terComponent;
	    terComponent.push_back(i);
	    terComponent.push_back(i + RequestNum);
	    sites.push_back(terComponent);
    }
}