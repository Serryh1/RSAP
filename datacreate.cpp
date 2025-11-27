#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <iomanip>
#include <string>

using namespace std;
struct Point
{
    double x, y;
};
int seed = 2333;
std::mt19937 gen(seed);

void generateData(int m, int k, int c, double sigma, int B, int groupId, int paramIndex, char paramChar)
{
    std::string filename_1 = "rasp_gaussian_m" + to_string(m) + "_k" + to_string(k) +
                             "_c" + to_string(c) + "_s" + to_string(static_cast<int>(sigma)) +
                             "_" + paramChar + to_string(paramIndex) + "_g" + to_string(groupId) + ".txt";

    string filename;
    if (paramChar == 'm')
    {
        filename = "../k_c_sigma/" + filename_1;
    }
    else if (paramChar == 'k')
    {
        filename = "../m_c_sigma/" + filename_1;
    }
    else if (paramChar == 'c')
    {
        filename = "../m_k_sigma/" + filename_1;
    }
    else if (paramChar == 's')
    {
        filename = "../m_k_c/" + filename_1;
    }

    int n = m * k;



    std::uniform_real_distribution<> uni_dist(0.0, B);           
    std::uniform_int_distribution<> center_dist(0, c - 1);       
    std::normal_distribution<> norm_dist(0.0, std::sqrt(sigma)); 


    std::vector<Point> centers;
    for (int i = 0; i < c; ++i)
    {
        centers.push_back({uni_dist(gen), uni_dist(gen)});
    }


    std::vector<Point> vehicles;
    for (int i = 0; i < m; ++i)
    {
        int j = center_dist(gen); 
        double x = centers[j].x + norm_dist(gen);
        double y = centers[j].y + norm_dist(gen);
        vehicles.push_back({x, y});
    }

    std::vector<Point> starts;
    for (int i = 0; i < n; ++i)
    {
        int j = center_dist(gen); 
        double x = centers[j].x + norm_dist(gen);
        double y = centers[j].y + norm_dist(gen);
        starts.push_back({x, y});
    }

    std::vector<Point> ends;
    for (int i = 0; i < n; ++i)
    {
        int j = center_dist(gen);
        double x = centers[j].x + norm_dist(gen);
        double y = centers[j].y + norm_dist(gen);
        ends.push_back({x, y});
    }

    std::ofstream out(filename);
    if (!out.is_open())
    {
        std::cerr << "Failed to open output file: " << filename << "\n";
        return;
    }

    out << "Vehicle " << m << "\n";
    for (auto &v : vehicles)
    {
        out << std::fixed << std::setprecision(2);
        out << "v " << v.x << " " << v.y << "\n";
    }

    out << "Request " << n << "\n";
    for (auto &s : starts)
    {
        out << std::fixed << std::setprecision(2);
        out << "s " << s.x << " " << s.y << "\n";
    }
    for (auto &t : ends)
    {
        out << std::fixed << std::setprecision(2);
        out << "t " << t.x << " " << t.y << "\n";
    }

    out.close();
    std::cout << "File written to " << filename << std::endl;
}

int main()
{
    int B = 4000; //  [0, B]

    vector<int> m_values = {100, 150, 200, 250, 300, 350, 400};     
    vector<int> k_values;         
    vector<int> c_values = {10, 15, 20, 25, 30, 35, 40};             
    vector<double> sigma_values = {40, 60, 80, 100, 120, 140, 160}; 

    for(int i = 2; i <= 34; ++ i){
        k_values.push_back(i);
    }


    int totalGroups = 10; 

    int fixed_m = m_values[3];
    int fixed_k = k_values[16];
    
    int fixed_c = c_values[3];
    double fixed_sigma = sigma_values[3];

    for (size_t m_idx = 0; m_idx < m_values.size(); ++m_idx)
    {
        int current_m = m_values[m_idx];
        for (int group = 1; group <= totalGroups; ++group)
        {
            generateData(current_m, fixed_k, fixed_c, fixed_sigma, B, group, m_idx + 1, 'm');
        }
    }


    for (size_t k_idx = 0; k_idx < k_values.size(); ++k_idx)
    {
        int current_k = k_values[k_idx];
        for (int group = 1; group <= totalGroups; ++group)
        {
            generateData(fixed_m, current_k, fixed_c, fixed_sigma, B, group, k_idx + 1, 'k');
        }
    }


    for (size_t c_idx = 0; c_idx < c_values.size(); ++c_idx)
    {
        int current_c = c_values[c_idx];
        for (int group = 1; group <= totalGroups; ++group)
        {
            generateData(fixed_m, fixed_k, current_c, fixed_sigma, B, group, c_idx + 1, 'c');
        }
    }

    for (size_t sigma_idx = 0; sigma_idx < sigma_values.size(); ++sigma_idx)
    {
        double current_sigma = sigma_values[sigma_idx];
        for (int group = 1; group <= totalGroups; ++group)
        {
            generateData(fixed_m, fixed_k, fixed_c, current_sigma, B, group, sigma_idx + 1, 's');
        }
    }

    cout << "All instances create success!" << endl;
    return 0;
}