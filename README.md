# RSAP

This is the experimental code of paper "Improved Algorithms for Trip-Vehicle Assignment in Ride-Sharing".

### Dataset

We have a total of four variables: m, k, c, $\sigma^{2}$. Each time, we vary one selected variable while keeping the other three fixed, sequentially generating datasets. The generated data is stored in the Dataset folder.

The Dataset folder contains four subfolders:

**k_c_sigma:** Here, k, c, and $\sigma^{2}$ remain constant, while m takes {100 150 200 **250** 300 350 400} 7 different values. 
**m_c_sigma:** Here, m, c, and $\sigma^{2}$ remain constant, while k takes {2  3 ... **18** ... 34} 33 different values.   
**m_k_c:** Here, m, k, and c remain constant, while $\sigma^{2}$  takes {40 60 80 **100** 120 140 160} 7 different values.  
**m_k_sigma:** Here, m, k, and $\sigma^{2}$ remain constant, while c takes {10 15 20 **25** 30 35 40} 7 different values.  

We fix three of the four parameters (m, k, c, $\sigma^{2}$) at their specified values, then generate results for all possible values of the remaining one.The data generation method is documented in the **datacreate.cpp**, and the random seed we use is **2333**.

The filenames follow the format:
	# rasp_gaussian_m100_k18_c25_s100_m1_g1.txt
where:

The first m (e.g., m100) indicates the value of m.

The second m (e.g., m1) denotes that this is the first value of m.

g1 signifies that this is the first dataset in this category.

For each fixed parameter combination, we generated 10 distinct datasets to ensure robustness.

The other folders use analogous naming conventions.

The format of the input file is as follows:  
1.  "Vehicle 100"： the number of vehicle;
2.  "v 3174.82 2381.07": the position of a vehicle;
3.  "Request 1800": the number of request;
4.  "s 2366.89 3228.28": the start position of the request;  
5.  "t 2339.47 3165.73": the end position of the request;  

### Environment

1. C++ 17;
2. cmake VERSION 3.5.0

### Run instructions
Create build folder

	mkdir build
Then use the following commands to generate the executable:

	cd build
	cmake .. 
	make 
Finally, execute the program using the following command, which requires three arguments:

	./program_name input_file_path output_file_path algorithm_name

an example is:

	cd ./build
	./RASP "../Dataset/k_c_sigma/rasp_gaussian_m100_k18_c25_s100_m1_g1.txt" "../output" "ALG1" 

### Expect output
This is an example output for instance rasp_gaussian_m250_k10_c25_s100_k9_g10.
#### 1. The output for ALG1 is as follows:
(1) after using the function GW_ALG1, it will output the running time of k-path construction.  
(2) after using the function Alg1_bigraphCreate and VehicleminCostMatching, it will output the group assigned to each vehicle. like:  

	Vehicle 181 choose group is : 0

(3) then using the function Alg1_computeResult, it will output the minumum cost path for each vehicle, like :  

	Vehicle number is :249
	minumum Cost is : 1776.365088
	The minumum cost Path is: 
	V249->S1158->S2438->S2044->S1711->S543->S2371->S1725->S6->S1572->S2165->T1158->T2438->T2044->T1711->T543->T2371->T1725->T6->T1572->T2165
V is the vehicle id, S and T is the start and end of the request.

(4) finally, it will output the total value and the running time, like:

	tot sum is : 851634.854160
	Running time is: 22.446868 s
	Computing the solution is complete!  

#### 2. The output for ALG2 is as follows:
(1) after using the function GW_ALG2, it will output the running time of k-path construction.  
(2) after using the function Alg2_bigraphCreate and VehicleminCostMatching, it will output the group assigned to each vehicle. like:  

	Vehicle 85 choose group is : 0
(3)  then using the function Alg2_computeResult, it will output the minimum cost path for each vehicle, like :

	Vehicle number is :249
	minumum Cost is : 4236.871140
	The minumum cost Path is: 
	V249->S486->S1724->S1854->S1572->S2461->S600->S1158->S2099->S2419->S543->T1158->T2419->T600->T2461->T1724->T2099->T543->T1854->T1572->T486

V is the vehicle id, S and T is the start and end of the request.
  
(4) finally, it will output the total value and the running time, like:

	TOT COST is: 1156422.434075
	Running time is: 99.911334 s
	Computing the solution is complete! 

#### 3. The output for HRA is as follow:
**Note on Algorithm Naming:**
While the algorithm is referred to as **LADG** in the paper, we implement it under the name **HRA** in our codebase. This nomenclature follows the original reference implementation by the algorithm's paper (who designated it as HRA in their source code).

(1) it will first output the pair of vehicle and group like follow:

	Vehicle 146 choose group is : 42
	Vehicle 166 choose group is : 43
	Vehicle 172 choose group is : 44
(2) then output each group's requset id:

	group id: 0 The requests are: 116 1866 127 1841 726 2005 1182 1524 0 576 
	group id: 1 The requests are: 197 1044 1298 2050 301 1483 302 1504 2 2128 
(3) then output each vehicle's path like follow:

	Start_S is :452
	End_S is :740
	Connect_T is :740
	....
	The final Path is: 
	V249->S452->S694->S760->S284->S561->S593->S363->S676->S698->S740->T740->T561->T694->T676->T593->T363->T698->T760->T284->T452

V is the vehicle id, S and T is the start and end of the request.

(4) finally, it will output the final result and the running time:

	TOT COST is: 1557628.352379
	Running time is: 22.575811 s
	Computing the solution is complete! 
