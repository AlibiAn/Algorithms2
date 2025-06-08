#ifndef TSP_SOLVER_H
#define TSP_SOLVER_H

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

// Core algorithm functions
std::vector<int> solve_tsp_cluster_hybrid(const std::vector<City>& cities, int K, int sub_problem_starts);
std::pair<std::vector<std::vector<int>>, std::vector<City>> kmeans_plus_plus(const std::vector<City>& cities, int K, int max_iter);
std::vector<int> solve_sub_problem(const std::vector<City>& sub_cities, int num_starts, const std::vector<City>& all_cities_ref);

// Connection and optimization functions
std::pair<int, int> find_closest_connection(const std::vector<int>& tour1, const std::vector<int>& tour2, const std::vector<City>& all_cities);
void apply_two_opt_local(std::vector<int>& local_tour, const std::vector<City>& sub_cities);
void apply_global_two_opt(std::vector<int>& tour, const std::vector<City>& all_cities);

// Tour construction heuristics
std::vector<int> nearest_neighbor_tour(const std::vector<City>& sub_cities, int start_node_local);
std::vector<int> furthest_insertion_start(const std::vector<City>& sub_cities);
std::vector<int> centroid_nearest_neighbor_start(const std::vector<City>& sub_cities);

// Utility functions
bool parse_tsplib_file(const std::string& filename, std::vector<City>& cities);
double distance(const City& c1, const City& c2);
double calculate_tour_cost(const std::vector<int>& tour, const std::vector<City>& all_cities);

#endif // TSP_SOLVER_H
