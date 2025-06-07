#include "result_tracker.h"
#include <iostream>
#include <iomanip>
#include <algorithm>

const std::map<std::string, double> ResultsTracker::KNOWN_OPTIMA = {
    {"a280", 2579.0},
    {"xql662", 2513.0},
    {"kz9976", 1061882.0},
    {"mona-lisa100K", 5757191.0}
};

ResultsTracker::ResultsTracker() {}

void ResultsTracker::addResult(const std::string& dataset, const std::string& algorithm, 
                              int vertices, double tour_cost, long long runtime_ms, 
                              const std::string& notes) {
    double optimal = getKnownOptimal(dataset);
    m_results.emplace_back(dataset, algorithm, vertices, tour_cost, runtime_ms, optimal, notes);
}

double ResultsTracker::getKnownOptimal(const std::string& dataset) {
    auto it = KNOWN_OPTIMA.find(dataset);
    return (it != KNOWN_OPTIMA.end()) ? it->second : -1.0;
}

void ResultsTracker::saveToCSV(const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot create file " << filename << std::endl;
        return;
    }
    
    file << "Dataset,Algorithm,Vertices,Tour_Cost,Runtime_ms,Known_Optimal,Approximation_Ratio,Notes\n";
    
    for (const auto& result : m_results) {
        file << result.m_dataset << ","
             << result.m_algorithm << ","
             << result.m_vertices << ","
             << std::fixed << std::setprecision(2) << result.m_tour_cost << ","
             << result.m_runtime_ms << ",";
        
        if (result.m_known_optimal > 0) {
            file << std::fixed << std::setprecision(2) << result.m_known_optimal << ","
                 << std::fixed << std::setprecision(3) << result.m_approximation_ratio;
        } else {
            file << "Unknown,N/A";
        }
        
        file << "," << result.m_notes << "\n";
    }
    
    file.close();
    std::cout << "Results saved to " << filename << std::endl;
}

void ResultsTracker::printSummary() {
    if (m_results.empty()) {
        std::cout << "No results to display" << std::endl;
        return;
    }
    
    std::cout << "\nTSP Algorithm Performance Summary" << std::endl;
    std::cout << std::setw(12) << "Dataset" 
              << std::setw(15) << "Algorithm" 
              << std::setw(8) << "Vertices"
              << std::setw(12) << "Cost"
              << std::setw(10) << "Time(ms)"
              << std::setw(8) << "Ratio" << std::endl;
    std::cout << std::string(65, '-') << std::endl;
    
    for (const auto& result : m_results) {
        std::cout << std::setw(12) << result.dataset
                  << std::setw(15) << result.algorithm
                  << std::setw(8) << result.vertices
                  << std::setw(12) << std::fixed << std::setprecision(0) << result.tour_cost
                  << std::setw(10) << result.runtime_ms;
        
        if (result.approximation_ratio > 0) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(2) << result.approximation_ratio;
        } else {
            std::cout << std::setw(8) << "N/A";
        }
        std::cout << std::endl;
    }
    std::cout << '\n' << std::endl;
}

void ResultsTracker::printComparison(const std::string& dataset) {
    std::cout << "\n Algorithm Comparison for " << dataset << std::endl;
    
    std::vector<TSPResult> dataset_results;
    for (const auto& result : m_results) {
        if (result.m_dataset == dataset) {
            dataset_results.push_back(result);
        }
    }
    
    if (dataset_results.empty()) {
        std::cout << "No results found for dataset: " << dataset << std::endl;
        return;
    }
    
    std::sort(dataset_results.begin(), dataset_results.end(), 
              [](const TSPResult& a, const TSPResult& b) {
                  return a.m_tour_cost < b.m_tour_cost;
              });
    
    std::cout << "Ranked by tour cost (best to worst):" << std::endl;
    for (size_t i = 0; i < dataset_results.size(); i++) {
        const auto& result = dataset_results[i];
        std::cout << (i + 1) << ". " << result.m_algorithm 
                  << " - Cost: " << std::fixed << std::setprecision(0) << result.m_tour_cost
                  << ", Time: " << result.m_runtime_ms << "ms";
        if (result.m_approximation_ratio > 0) {
            std::cout << ", Ratio: " << std::fixed << std::setprecision(2) << result.m_approximation_ratio;
        }
        std::cout << std::endl;
    }
    std::cout << '\n' << std::endl;
}

void ResultsTracker::clear() {
    m_results.clear();
}

std::vector<TSPResult> ResultsTracker::getResults() const {
    return m_results;
}
