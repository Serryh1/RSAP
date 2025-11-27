#include "SF-Solution.h"
#include "SF-Functions.h"
#include "GW.h"
#include "Input.h"
#include "BipartiteGraph.h"
#include "HRA.h"

#include <time.h>
#include <sys/time.h>
#include <sstream>
#include <string>
#include <string.h>
#include <fstream>
#include <cassert>
#include <chrono>


using namespace std;

int main(int argc, char* argv[])
{
    
    string alg,strInputFile,strOutputFile;

    int nVertex,nEdges,nTerminals;
    vector<vector<double>> adjMatrix;
    vector<vector<int>> E;
    int *terminals;
    vector<vector<int> > pairedTerminals;
    
    int VehicleNum;
    int RequestNum;
    int K;

    vector<pair<double,double>> Rs;
    vector<pair<double,double>> Rt;
    vector<pair<double,double>> VehilcePosition;
    vector<vector<int>> kPath;

    if(argc < 3){
        cout << "Not enough parameters have been passed. \n";
        cin.get();
        exit(0);
    }
    else{
        strInputFile=argv[1];
        strOutputFile=argv[2];
        alg = argv[3];
    }

    // strInputFile = "../rasp_gaussian_m250_k6_c25_s100_k5_g6.txt";
    // strOutputFile = "output.txt";
    // alg = "HRA";

    std::ofstream file(strOutputFile);
    std::cout.rdbuf(file.rdbuf());  

	cout<<"strInputFile: "<<strInputFile<<endl;
	cout<<"strOutputFile: "<<strOutputFile<<endl;

    if (alg == "ALG1"){

        cout << "The Algorithm is for ALG1" << endl;
        cout<<"********************* File  "<< strInputFile <<"  is read for GW algorithms. ********************** \n";
        cout<<"Calling GW_ALG1 Algorithm \n";

        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();

        ALG1_InputfileProcess(strInputFile, nVertex, nEdges, nTerminals, terminals, adjMatrix, pairedTerminals, E, VehicleNum, RequestNum, Rs, Rt, VehilcePosition);
        
        assert(RequestNum % VehicleNum == 0);
        K = RequestNum / VehicleNum;
        cout << K << '\n';

        GW_ALG1(nVertex, nEdges, nTerminals, terminals, adjMatrix, pairedTerminals, E, K, kPath);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "The running time for GW and split is : " << duration.count() << " s" << std::endl;

        // BipartiteGraph Create
        vector<vector<double>> BipartiteGraphAdj; 
        vector<int> BipartiteGraphResult; // vehicle ---> groupid
        Alg1_bigraphCreate(VehilcePosition, kPath, Rs, Rt, BipartiteGraphAdj);
        VehicleminCostMatching(BipartiteGraphAdj, BipartiteGraphResult);
        Alg1_computeResult(VehilcePosition, kPath, Rs, Rt, BipartiteGraphResult);
        end = std::chrono::high_resolution_clock::now();
        duration = end - start;
        std::cout << "Running time is: " << duration.count() << " s" << std::endl;

    }else if(alg == "ALG2"){

        cout << "The Algorithm is for ALG2" << endl;
        cout<<"********************* File  "<< strInputFile <<"  is read for GW algorithms. ********************** \n";
        cout<<"Calling GW_ALG2 Algorithm \n";

        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();
        ALG2_InputfileProcess(strInputFile, nVertex, nEdges, nTerminals, terminals, adjMatrix, pairedTerminals, E, VehicleNum, RequestNum, Rs, Rt, VehilcePosition);
        
        // look for 2k path
        assert(RequestNum % VehicleNum == 0);
        K = RequestNum / VehicleNum * 2;
        GW_ALG2(nVertex, nEdges, nTerminals, terminals, adjMatrix, pairedTerminals, E, K, kPath);
        // BipartiteGraph Create
        vector<vector<double>> BipartiteGraphAdj; 
        vector<int> BipartiteGraphResult; // vehicle ---> groupid
        Alg2_bigraphCreate(VehilcePosition, kPath, Rs, Rt, BipartiteGraphAdj, nVertex / 2);
        VehicleminCostMatching(BipartiteGraphAdj, BipartiteGraphResult);
        Alg2_computeResult(VehilcePosition, kPath, Rs, Rt, BipartiteGraphResult, RequestNum);

        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "Running time is: " << duration.count() << " s" << std::endl;
    }else if(alg == "HRA"){

        cout << "The Algorithm is for HRA" << endl;
        cout<<"********************* File  "<< strInputFile <<"  is read for HRA algorithms. ********************** \n";
        cout<<"Calling HRA Algorithm \n";

        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();
        ALG1_InputfileProcess(strInputFile, nVertex, nEdges, nTerminals, terminals, adjMatrix, pairedTerminals, E, VehicleNum, RequestNum, Rs, Rt, VehilcePosition);
        assert(RequestNum % VehicleNum == 0);
        
        K = RequestNum / VehicleNum;
        HRA Hra;
        Hra.HRASolve(K, VehilcePosition, Rs, Rt, RequestNum, VehicleNum);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "Running time is: " << duration.count() << " s" << std::endl;

    }else{
        cout << "No such an algorithm!" << endl;
        return 0;
    }
    
    cout<<"Computing the solution is complete! \n";
    cout<<"--------------------------------------------"<<endl;
    cout << strInputFile << "\n" << nVertex <<"\n"<< nEdges << "\n" << nTerminals <<"\n";
    cout<<"******************************************************************** \n";
    return 0;
}

