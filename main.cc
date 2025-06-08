#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <map>

#include "mst.h"
#include "karp.h"
#include "mine.h"

//struct to hold results for CSV output
struct TSPResult {
    std::string dataset;
    std::string algorithm;
    int vertices;
    double cost;
    long long runtime_ms;
    double approx_ratio;
};

// Optimal costs from the website()
const std::map<std::string, double> KNOWN_OPTIMA = {
    {"a280", 2579.0}, {"xql662", 2513.0}, {"kz9976", 1061882.0}, {"mona-lisa100K", 5757191.0}
};

void saveResultsToCSV(const std::vector<TSPResult>& results, const std::string& filename) {
    std::ofstream file(filename);
    file << "Dataset,Algorithm,Vertices,Tour_Cost,Runtime_ms,Approximation_Ratio\n";
    for (const auto& res : results) {
        file << res.dataset << "," << res.algorithm << "," << res.vertices << ","
             << std::fixed << std::setprecision(2) << res.cost << "," << res.runtime_ms << ",";
        if (res.approx_ratio > 0) {
            file << std::fixed << std::setprecision(4) << res.approx_ratio;
        }
        file << "\n";
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "proper usage: " << argv[0] << " <dataset1.tsp> [dataset2.tsp] ..." << std::endl;
        return 1;
    }

    std::vector<TSPResult> all_results;
    std::vector<std::string> datasets(argv + 1, argv + argc);

    for (const auto& filename : datasets) {
        std::cout << "\n Processing Dataset: " << filename << std::endl;
        
        //runs MstTsp
        MSTTSPSolver mst_solver;
        if (mst_solver.loadFromFile(filename)) {
            auto start = std::chrono::high_resolution_clock::now();
            auto [tour, cost] = mst_solver.solve();
            auto end = std::chrono::high_resolution_clock::now();
            long long duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            
            double ratio = -1.0;
            if (KNOWN_OPTIMA.count(mst_solver.getInstanceName())) {
                ratio = cost / KNOWN_OPTIMA.at(mst_solver.getInstanceName());
            }
            all_results.push_back({mst_solver.getInstanceName(), "MST-2-Approx", mst_solver.getVertexCount(), cost, duration, ratio});
            std::cout << " MST Cost: " << cost << ", Time: " << duration << "ms" << std::endl;
        }

        //RUNs held-karp
        HeldKarpTSPSolver karp_solver;
        if (karp_solver.loadFromFile(filename) && karp_solver.getVertexCount() <= 25) {
            auto start = std::chrono::high_resolution_clock::now();
            auto [tour, cost] = karp_solver.solve();
            auto end = std::chrono::high_resolution_clock::now();
            long long duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            
            double ratio = -1.0;
            if (cost > 0 && KNOWN_OPTIMA.count(karp_solver.getInstanceName())) {
                ratio = cost / KNOWN_OPTIMA.at(karp_solver.getInstanceName());
            }
            all_results.push_back({karp_solver.getInstanceName(), "Held-Karp", karp_solver.getVertexCount(), cost, duration, ratio});
            std::cout << "Held-karp cost: " << cost << ", Time: " << duration << "ms" << std::endl;
        } else if(karp_solver.getVertexCount() > 25) {
             std::cout << "Skipping Held-Karp: Too many vertices for my Hardware (" << karp_solver.getVertexCount() << ")." << std::endl;
        }

        //my wonderful algorithm
        std::vector<City> cities;
        if (parse_tsplib_file(filename, cities)) {
            auto start = std::chrono::high_resolution_clock::now();
            std::vector<int> tour = solve_tsp_cluster_hybrid(cities, 20, 10);
            double cost = calculate_tour_cost(tour, cities);
            auto end = std::chrono::high_resolution_clock::now();
            long long duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

            double ratio = -1.0;
            std::string dataset_name = filename.substr(filename.find_last_of("/\\") + 1);
            dataset_name = dataset_name.substr(0, dataset_name.find(".tsp"));
            if (KNOWN_OPTIMA.count(dataset_name)) {
                ratio = cost / KNOWN_OPTIMA.at(dataset_name);
            }
            all_results.push_back({dataset_name, "Mine", (int)cities.size(), cost, duration, ratio});
            std::cout << "My algorithm's cost: " << cost << ", Time: " << duration << "ms" << std::endl;
        }
    }

    saveResultsToCSV(all_results, "tsp_results.csv");

    return 0;
}
