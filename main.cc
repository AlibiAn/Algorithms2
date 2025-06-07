#include "mst.h"
#include "held_karp.h"
#include "results_tracker.h"
#include <iostream>
#include <vector>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <dataset_file> [--csv]" << std::endl;
        return 1;
    }
    
    std::string filename = argv[1];
    bool csv_output = (argc > 2 && std::string(argv[2]) == "--csv");
    
    ResultsTracker tracker;
    
    // Extract dataset name from filename
    std::string dataset = filename.substr(filename.find_last_of("/\\") + 1);
    if (dataset.find(".tsp") != std::string::npos) {
        dataset = dataset.substr(0, dataset.find(".tsp"));
    }
    
    // Test MST Algorithm
    MSTTSPSolver mst_solver;
    if (mst_solver.loadFromFile(filename)) {
        auto [mst_tour, mst_cost, mst_time] = mst_solver.solveWithTiming();
        tracker.addResult(dataset, "MST-2-Approx", mst_solver.getVertexCount(), 
                         mst_cost, mst_time, "2-approximation");
        
        if (!csv_output) {
            mst_solver.printTourResults(mst_tour, mst_cost, mst_time);
        }
    }
    
    // Test Held-Karp Algorithm (only for small instances)
    if (mst_solver.getVertexCount() <= 25) {
        HeldKarpTSPSolver hk_solver;
        if (hk_solver.loadFromFile(filename)) {
            auto [hk_tour, hk_cost, hk_time] = hk_solver.solveWithTiming();
            tracker.addResult(dataset, "Held-Karp", hk_solver.getVertexCount(), 
                             hk_cost, hk_time, "Exact optimal");
            
            if (!csv_output) {
                hk_solver.printTourResults(hk_tour, hk_cost, hk_time);
            }
        }
    } else {
        if (!csv_output) {
            std::cout << "Skipping Held-Karp for large instance (" 
                      << mst_solver.getVertexCount() << " vertices)" << std::endl;
        }
    }
    
    // Add your novel algorithm here when implemented
    // tracker.addResult(dataset, "My-Algorithm", vertices, cost, time, "Novel approach");
    
    if (csv_output) {
        // Output in CSV format for the Makefile
        for (const auto& result : tracker.getResults()) {
            std::cout << result.dataset << "," << result.algorithm << "," 
                      << result.vertices << "," << result.tour_cost << "," 
                      << result.runtime_ms << ",";
            if (result.approximation_ratio > 0) {
                std::cout << result.approximation_ratio;
            } else {
                std::cout << "N/A";
            }
            std::cout << std::endl;
        }
    } else {
        tracker.printSummary();
        tracker.printComparison(dataset);
        tracker.saveToCSV("results_" + dataset + ".csv");
    }
    
    return 0;
}
