#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <numeric>
#include <cmath>

#include "mst.h"
#include "karp.h"
#include "mine.h"

struct TSPResult {
    std::string dataset;
    std::string algorithm;
    int vertices;
    double cost;
    long long runtime_ms;
    double approx_ratio;
    int k_clusters;
    int num_starts;
    int run_id;
};

const std::map<std::string, double> KNOWN_OPTIMA = {
    {"a280", 2579.0}, 
    {"xql662", 2513.0}, 
    {"kz9976", 1061882.0}, 
    {"mona-lisa100K", 5757191.0}
};

int getDatasetSize(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return -1;
    
    std::string line;
    while (std::getline(file, line)) {
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (line.empty()) continue;
        
        if (line.find("DIMENSION") != std::string::npos) {
            size_t colon_pos = line.find(':');
            std::string dim_str = (colon_pos != std::string::npos) ? line.substr(colon_pos + 1) : line.substr(9);
            dim_str.erase(0, dim_str.find_first_not_of(" \t"));
            try { return std::stoi(dim_str); } catch (...) { return -1; }
        }
    }
    return -1;
}

void saveResultsToCSV(const std::vector<TSPResult>& results, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot create " << filename << std::endl;
        return;
    }
    
    file << "Dataset,Algorithm,Vertices,K_Clusters,Num_Starts,Run_ID,Tour_Cost,Runtime_ms,Approximation_Ratio\n";
    for (const auto& res : results) {
        file << res.dataset << "," << res.algorithm << "," << res.vertices << ","
             << res.k_clusters << "," << res.num_starts << "," << res.run_id << ","
             << std::fixed << std::setprecision(2) << res.cost << "," << res.runtime_ms << ",";
        if (res.approx_ratio > 0) file << std::fixed << std::setprecision(4) << res.approx_ratio;
        file << "\n";
    }
    std::cout << "\nAnalysis results saved to " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <dataset1.tsp> [dataset2.tsp] ..." << std::endl;
        return 1;
    }

    std::vector<TSPResult> all_results;
    std::vector<std::string> datasets(argv + 1, argv + argc);

    const int MST_KARP_RUNS = 10;
    const int MINE_RUNS_PER_PARAM = 5;

    for (const auto& filename : datasets) {
        std::cout << "\n--- Processing Dataset: " << filename << " ---" << std::endl;
        
        int dataset_size = getDatasetSize(filename);
        if (dataset_size == -1) {
            std::cout << "Skipping " << filename << ": Cannot read dataset size." << std::endl;
            continue;
        }
        std::cout << "Dataset size: " << dataset_size << " vertices" << std::endl;
        
        std::string dataset_name = filename.substr(filename.find_last_of("/\\") + 1);
        if (dataset_name.find(".tsp") != std::string::npos) {
            dataset_name = dataset_name.substr(0, dataset_name.find(".tsp"));
        }

        std::cout << "Testing MST Algorithm (" << MST_KARP_RUNS << " runs)..." << std::endl;
        for (int run = 1; run <= MST_KARP_RUNS; ++run) {
            MSTTSPSolver mst_solver;
            if (mst_solver.loadFromFile(filename)) {
                auto start_time = std::chrono::high_resolution_clock::now();
                auto [tour, cost] = mst_solver.solve();
                auto end_time = std::chrono::high_resolution_clock::now();
                long long duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
                double ratio = KNOWN_OPTIMA.count(dataset_name) ? cost / KNOWN_OPTIMA.at(dataset_name) : -1.0;
                all_results.push_back({dataset_name, "MST-2-Approx", dataset_size, cost, duration, ratio, -1, -1, run});
            }
        }

        if (dataset_size <= 25) {
            std::cout << "Testing Held-Karp Algorithm (" << MST_KARP_RUNS << " runs)..." << std::endl;
            for (int run = 1; run <= MST_KARP_RUNS; ++run) {
                HeldKarpTSPSolver karp_solver;
                if (karp_solver.loadFromFile(filename)) {
                    auto start_time = std::chrono::high_resolution_clock::now();
                    auto [tour, cost] = karp_solver.solve();
                    auto end_time = std::chrono::high_resolution_clock::now();
                    long long duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
                    double ratio = (cost > 0 && KNOWN_OPTIMA.count(dataset_name)) ? cost / KNOWN_OPTIMA.at(dataset_name) : -1.0;
                    all_results.push_back({dataset_name, "Held-Karp", dataset_size, cost, duration, ratio, -1, -1, run});
                }
            }
        }

        if (dataset_size > 15000) {
            std::cout << "Skipping My Algorithm: Dataset too large for detailed analysis." << std::endl;
        } else {
            std::vector<City> cities;
            if (parse_tsplib_file(filename, cities)) {
                int k_mid = static_cast<int>(sqrt(dataset_size));
                std::vector<int> k_values = {std::max(2, k_mid / 2), k_mid, k_mid * 2};

                for (int k : k_values) {
                    std::cout << "Testing My Algorithm with K=" << k << " (" << MINE_RUNS_PER_PARAM << " runs)..." << std::endl;
                    for (int run = 1; run <= MINE_RUNS_PER_PARAM; ++run) {
                        auto start_time = std::chrono::high_resolution_clock::now();
                        std::vector<int> tour = solve_tsp_cluster_hybrid(cities, k, 10);
                        double cost = calculate_tour_cost(tour, cities);
                        auto end_time = std::chrono::high_resolution_clock::now();
                        long long duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
                        double ratio = KNOWN_OPTIMA.count(dataset_name) ? cost / KNOWN_OPTIMA.at(dataset_name) : -1.0;
                        all_results.push_back({dataset_name, "Mine", dataset_size, cost, duration, ratio, k, 10, run});
                    }
                }
            }
        }
    }

    saveResultsToCSV(all_results, "tsp_results.csv");
    return 0;
}
