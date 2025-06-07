#ifndef RESULTS_TRACKER_H
#define RESULTS_TRACKER_H

#include <string>
#include <vector>
#include <fstream>
#include <map>

struct TSPResult {
    std::string m_dataset;
    std::string m_algorithm;
    int m_vertices;
    double m_tour_cost;
    long long m_runtime_ms;
    double m_known_optimal;
    double m_approximation_ratio;
    std::string m_notes;
    
    TSPResult(const std::string& ds, const std::string& alg, int n, 
              double cost, long long time, double optimal = -1.0, const std::string& note = "")
        : m_dataset(ds), m_algorithm(alg), m_vertices(n), m_tour_cost(cost), 
          m_runtime_ms(time), m_known_optimal(optimal), m_notes(note) {
        if (optimal > 0) {
            m_approximation_ratio = cost / optimal;
        } else {
            m_approximation_ratio = -1.0;
        }
    }
};

class ResultsTracker {
private:
    std::vector<TSPResult> m_results;
    static const std::map<std::string, double> KNOWN_OPTIMA;
    
public:
    ResultsTracker();
    
    void addResult(const std::string& dataset, const std::string& algorithm, 
                   int vertices, double tour_cost, long long runtime_ms, 
                   const std::string& notes = "");
    
    void saveToCSV(const std::string& filename = "tsp_results.csv");
    void printSummary();
    void printComparison(const std::string& dataset);
    void clear();
    
    double getKnownOptimal(const std::string& dataset);
    std::vector<TSPResult> getResults() const;
};

#endif
