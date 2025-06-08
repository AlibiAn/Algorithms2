#ifndef MINE_H
#define MINE_H

#include <vector>
#include <string>
#include <utility>

//// Version 3.1: Some Performance Improvements!
// Changes from V3:
//1)Enhanced 2-opt with epsilon and early termination
//2)Post-processing 2-opt on final tour  
//3)Mixed starting strategies for sub-problems
//4)Adaptive effort based on cluster size
//5)Bidirectional stitching evaluation
//
// Notes: Slightly better than V3.0 but does not scale very well too

struct City { 
    int id; 
    double x, y; 
};

std::vector<int> solve_tsp_cluster_hybrid(const std::vector<City>& cities, int K, int sub_problem_starts);

bool parse_tsplib_file(const std::string& filename, std::vector<City>& cities);
double calculate_tour_cost(const std::vector<int>& tour, const std::vector<City>& all_cities);

#endif
