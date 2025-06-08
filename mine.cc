#include "mine.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <numeric>
#include <limits>
#include <algorithm>
#include <random>
#include <map>

std::vector<int> solve_tsp_cluster_hybrid(const std::vector<City>& cities, int K, int sub_problem_starts) {
    std::cout << "\nclustering cities" << std::endl;
    auto [clusters_of_global_ids, centroids] = kmeans_pp(cities, K, 100);

    std::cout << "\n solving TSP for each cluster" << std::endl;
    std::vector<std::vector<int>> sub_tours_global(K);
    for (int i = 0; i < K; ++i) {
        if (clusters_of_global_ids[i].empty()) continue;
        std::vector<City> sub_problem_cities;
        for (int global_id : clusters_of_global_ids[i]) {
            sub_problem_cities.push_back(cities[global_id]);
        }
        int adaptive_starts = std::max(3, std::min(sub_problem_starts, static_cast<int>(sub_problem_cities.size() * 1.5)));
        sub_tours_global[i] = solve_sub_problem(sub_problem_cities, adaptive_starts, cities);
    }
    
    std::cout << "\n finding the optimal cluster sequence" << std::endl;
    std::vector<int> centroid_tour_local = solve_sub_problem(centroids, K * 2, centroids);
    if (centroid_tour_local.empty()) return {};
    centroid_tour_local.pop_back();

    std::cout << "\n stitching together the sub-tours with the closest node pair" << std::endl;
    if (centroid_tour_local.empty()) return {};

    int first_cluster_idx = centroid_tour_local[0];
    std::vector<int> final_tour = sub_tours_global[first_cluster_idx];
    if (final_tour.empty()) return {};
    final_tour.pop_back();

    for (size_t i = 1; i < centroid_tour_local.size(); ++i) {
        int next_cluster_idx = centroid_tour_local[i];
        std::vector<int> next_sub_tour = sub_tours_global[next_cluster_idx];
        if (next_sub_tour.empty()) continue;
        next_sub_tour.pop_back();

        auto [node_from_final, node_to_attach] = find_closest_connection(final_tour, next_sub_tour, cities);

        auto it_final = std::find(final_tour.begin(), final_tour.end(), node_from_final);
        if (it_final != final_tour.end()) {
             std::rotate(final_tour.begin(), it_final + 1, final_tour.end());
        }
        
        auto it_next = std::find(next_sub_tour.begin(), next_sub_tour.end(), node_to_attach);
        if (it_next != next_sub_tour.end()) {
            std::rotate(next_sub_tour.begin(), it_next, next_sub_tour.end());
        }
        
        final_tour.insert(final_tour.end(), next_sub_tour.begin(), next_sub_tour.end());
    }

    if (!final_tour.empty()) final_tour.push_back(final_tour.front());

    std::cout << "\nnew: post-processing with global 2-opt" << std::endl;
    apply_global_two_opt(final_tour, cities);

    return final_tour;
}

std::pair<int, int> find_closest_connection(const std::vector<int>& tour1, const std::vector<int>& tour2, const std::vector<City>& all_cities) {
    int best_node1 = -1, best_node2 = -1;
    double min_dist = std::numeric_limits<double>::max();

    if (tour1.empty() || tour2.empty()) return {-1, -1};

    // Try normal orientation
    for (int node1_id : tour1) {
        for (int node2_id : tour2) {
            double d = distance(all_cities[node1_id], all_cities[node2_id]);
            if (d < min_dist) {
                min_dist = d;
                best_node1 = node1_id;
                best_node2 = node2_id;
            }
        }
    }

    // Try reversed orientation of tour2
    std::vector<int> tour2_rev(tour2.rbegin(), tour2.rend());
    for (int node1_id : tour1) {
        for (int node2_id : tour2_rev) {
            double d = distance(all_cities[node1_id], all_cities[node2_id]);
            if (d < min_dist) {
                min_dist = d;
                best_node1 = node1_id;
                best_node2 = node2_id;
            }
        }
    }
    
    return {best_node1, best_node2};
}

void apply_two_opt_local(std::vector<int>& local_tour, const std::vector<City>& sub_cities) {
    bool improved = true;
    int n = static_cast<int>(local_tour.size()) - 1;
    const double epsilon = 1e-9;
    
    while (improved) {
        improved = false;
        for (int i = 0; i < n - 1; ++i) {
            for (int j = i + 1; j < n; ++j) {
                const City& c_i = sub_cities[local_tour[i]];
                const City& c_i1 = sub_cities[local_tour[i+1]];
                const City& c_j = sub_cities[local_tour[j]];
                const City& c_j1 = sub_cities[local_tour[j+1]];
                
                double old_dist = distance(c_i, c_i1) + distance(c_j, c_j1);
                double new_dist = distance(c_i, c_j) + distance(c_i1, c_j1);
                
                if (new_dist < old_dist - epsilon) {
                    std::reverse(local_tour.begin() + i + 1, local_tour.begin() + j + 1);
                    improved = true;
                    break;
                }
            }
            if (improved) break;
        }
    }
}

void apply_global_two_opt(std::vector<int>& tour, const std::vector<City>& all_cities) {
    bool improved = true;
    int n = static_cast<int>(tour.size()) - 1;
    const double epsilon = 1e-9;
    int iterations = 0;
    const int max_iterations = 3;
    
    while (improved && iterations < max_iterations) {
        improved = false;
        iterations++;
        
        for (int i = 0; i < n - 1; ++i) {
            for (int j = i + 1; j < n; ++j) {
                const City& c_i = all_cities[tour[i]];
                const City& c_i1 = all_cities[tour[i+1]];
                const City& c_j = all_cities[tour[j]];
                const City& c_j1 = all_cities[tour[j+1]];
                
                double old_dist = distance(c_i, c_i1) + distance(c_j, c_j1);
                double new_dist = distance(c_i, c_j) + distance(c_i1, c_j1);
                
                if (new_dist < old_dist - epsilon) {
                    std::reverse(tour.begin() + i + 1, tour.begin() + j + 1);
                    improved = true;
                    break;
                }
            }
            if (improved) break;
        }
    }
}

std::vector<int> solve_sub_problem(const std::vector<City>& sub_cities, int num_starts, const std::vector<City>& all_cities_ref) {
    const int num_sub_cities = static_cast<int>(sub_cities.size());
    if (num_sub_cities == 0) return {};
    if (num_sub_cities == 1) return {sub_cities[0].id, sub_cities[0].id};

    std::vector<int> best_tour_local;
    double min_cost = std::numeric_limits<double>::max();
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0, num_sub_cities - 1);

    for (int i = 0; i < num_starts; ++i) {
        std::vector<int> current_tour_local;
        
        if (i == 0 && num_sub_cities > 2) {
            current_tour_local = furthest_insertion_start(sub_cities);
        } else if (i == 1 && num_sub_cities > 3) {
            current_tour_local = centroid_nearest_neighbor_start(sub_cities);
        } else {
            int start_node_local = distrib(gen);
            current_tour_local = nearest_neighbor_tour(sub_cities, start_node_local);
        }
        
        if (current_tour_local.empty()) continue;
        
        apply_two_opt_local(current_tour_local, sub_cities);

        std::vector<int> temp_global_tour;
        for(int local_idx : current_tour_local) temp_global_tour.push_back(sub_cities[local_idx].id);
        
        double current_cost = calculate_tour_cost(temp_global_tour, all_cities_ref);
        if (current_cost < min_cost) {
            min_cost = current_cost;
            best_tour_local = current_tour_local;
        }
    }
    
    std::vector<int> best_tour_global;
    for(int local_idx : best_tour_local) best_tour_global.push_back(sub_cities[local_idx].id);
    return best_tour_global;
}

std::vector<int> nearest_neighbor_tour(const std::vector<City>& sub_cities, int start_node_local) {
    const int num_sub_cities = static_cast<int>(sub_cities.size());
    std::vector<int> current_tour_local;
    std::vector<bool> visited(num_sub_cities, false);
    current_tour_local.push_back(start_node_local);
    visited[start_node_local] = true;

    while(static_cast<int>(current_tour_local.size()) < num_sub_cities){
        int best_next_local = -1;
        double min_dist = std::numeric_limits<double>::max();
        for(int next_local = 0; next_local < num_sub_cities; ++next_local){
            if(!visited[next_local]){
                double d = distance(sub_cities[current_tour_local.back()], sub_cities[next_local]);
                if(d < min_dist){
                    min_dist = d;
                    best_next_local = next_local;
                }
            }
        }
        if(best_next_local != -1) { 
            current_tour_local.push_back(best_next_local); 
            visited[best_next_local] = true; 
        } else break;
    }
    current_tour_local.push_back(start_node_local);
    return current_tour_local;
}

std::vector<int> furthest_insertion_start(const std::vector<City>& sub_cities) {
    const int n = static_cast<int>(sub_cities.size());
    if (n < 3) return nearest_neighbor_tour(sub_cities, 0);
    
    double max_dist = 0;
    int city1 = 0, city2 = 1;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double d = distance(sub_cities[i], sub_cities[j]);
            if (d > max_dist) {
                max_dist = d;
                city1 = i;
                city2 = j;
            }
        }
    }
    
    std::vector<int> tour = {city1, city2, city1};
    std::vector<bool> in_tour(n, false);
    in_tour[city1] = in_tour[city2] = true;
    
    while (static_cast<int>(tour.size()) - 1 < n) {
        int best_city = -1;
        double max_min_dist = -1;
        
        for (int i = 0; i < n; ++i) {
            if (in_tour[i]) continue;
            double min_dist_to_tour = std::numeric_limits<double>::max();
            for (int j = 0; j < static_cast<int>(tour.size()) - 1; ++j) {
                min_dist_to_tour = std::min(min_dist_to_tour, distance(sub_cities[i], sub_cities[tour[j]]));
            }
            if (min_dist_to_tour > max_min_dist) {
                max_min_dist = min_dist_to_tour;
                best_city = i;
            }
        }
        
        if (best_city == -1) break;
        
        int best_pos = 1;
        double min_increase = std::numeric_limits<double>::max();
        for (int pos = 1; pos < static_cast<int>(tour.size()); ++pos) {
            double increase = distance(sub_cities[tour[pos-1]], sub_cities[best_city]) + 
                            distance(sub_cities[best_city], sub_cities[tour[pos]]) -
                            distance(sub_cities[tour[pos-1]], sub_cities[tour[pos]]);
            if (increase < min_increase) {
                min_increase = increase;
                best_pos = pos;
            }
        }
        
        tour.insert(tour.begin() + best_pos, best_city);
        in_tour[best_city] = true;
    }
    
    return tour;
}

std::vector<int> centroid_nearest_neighbor_start(const std::vector<City>& sub_cities) {
    const int n = static_cast<int>(sub_cities.size());
    
    double cx = 0, cy = 0;
    for (const auto& city : sub_cities) {
        cx += city.x;
        cy += city.y;
    }
    cx /= n;
    cy /= n;
    
    int closest_to_centroid = 0;
    double min_dist = std::numeric_limits<double>::max();
    for (int i = 0; i < n; ++i) {
        double d = std::sqrt((sub_cities[i].x - cx) * (sub_cities[i].x - cx) + 
                           (sub_cities[i].y - cy) * (sub_cities[i].y - cy));
        if (d < min_dist) {
            min_dist = d;
            closest_to_centroid = i;
        }
    }
    
    return nearest_neighbor_tour(sub_cities, closest_to_centroid);
}

std::pair<std::vector<std::vector<int>>, std::vector<City>> kmeans_pp(const std::vector<City>& cities, int K, int max_iter) {
    const int N = static_cast<int>(cities.size());
    if (K > N) K = N;
    std::vector<City> centroids;
    std::random_device rd;
    std::mt19937 gen(rd());
    
    std::uniform_int_distribution<> distrib(0, N - 1);
    centroids.push_back(cities[distrib(gen)]);
    
    std::vector<double> dist_sq(N);
    for (int k = 1; k < K; ++k) {
        double total_dist_sq = 0;
        for (int i = 0; i < N; ++i) {
            double min_dist_sq = std::numeric_limits<double>::max();
            for (const auto& centroid : centroids) {
                double d = distance(cities[i], centroid);
                min_dist_sq = std::min(min_dist_sq, d * d);
            }
            dist_sq[i] = min_dist_sq;
            total_dist_sq += dist_sq[i];
        }
        if(total_dist_sq == 0) continue;
        std::uniform_real_distribution<> r_dist(0, total_dist_sq);
        double r = r_dist(gen);
        double cumulative_dist = 0;
        for (int i = 0; i < N; ++i) {
            cumulative_dist += dist_sq[i];
            if (cumulative_dist >= r) {
                centroids.push_back(cities[i]);
                break;
            }
        }
    }

    std::vector<int> assignments(N);
    for (int iter = 0; iter < max_iter; ++iter) {
        for (int i = 0; i < N; ++i) {
            double min_dist = std::numeric_limits<double>::max();
            int best_cluster = 0;
            for (int k = 0; k < static_cast<int>(centroids.size()); ++k) {
                double d = distance(cities[i], centroids[k]);
                if (d < min_dist) { min_dist = d; best_cluster = k; }
            }
            assignments[i] = best_cluster;
        }
        std::vector<City> new_centroids(K, {0, 0.0, 0.0});
        std::vector<int> counts(K, 0);
        for (int i = 0; i < N; ++i) {
            new_centroids[assignments[i]].x += cities[i].x;
            new_centroids[assignments[i]].y += cities[i].y;
            counts[assignments[i]]++;
        }
        bool converged = true;
        for (int k = 0; k < K; ++k) {
            if (counts[k] > 0) {
                new_centroids[k].x /= counts[k];
                new_centroids[k].y /= counts[k];
                if (distance(new_centroids[k], centroids[k]) > 1e-6) converged = false;
                centroids[k] = new_centroids[k];
            }
        }
        if (converged) break;
    }

    std::vector<std::vector<int>> clusters(K);
    for (int i = 0; i < N; ++i) clusters[assignments[i]].push_back(cities[i].id);
    for(int k=0; k<K; ++k) centroids[k].id = k;
    return {clusters, centroids};
}

bool parse_tsplib_file(const std::string& filename, std::vector<City>& cities) {
    std::ifstream file(filename);
    if (!file.is_open()) { 
        std::cerr << "cerrror: Could'nt open file " << filename << std::endl; 
        return false; 
    }
    std::string line;
    bool in_coord_section = false;
    while (std::getline(file, line)) {
        if (line.find("NODE_COORD_SECTION") != std::string::npos) { 
            in_coord_section = true; 
            continue; 
        }
        if (line.find("EOF") != std::string::npos) break;
        if (in_coord_section) {
            std::stringstream ss(line);
            City city;
            if (ss >> city.id >> city.x >> city.y) {
                city.id--;
                cities.push_back(city);
            }
        }
    }
    std::sort(cities.begin(), cities.end(), [](const City& a, const City& b){ return a.id < b.id; });
    return true;
}

double distance(const City& c1, const City& c2) {
    return std::round(std::sqrt(std::pow(c1.x - c2.x, 2) + std::pow(c1.y - c2.y, 2)));
}

double calculate_tour_cost(const std::vector<int>& tour, const std::vector<City>& all_cities) {
    if (tour.empty()) return 0.0;
    double cost = 0.0;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
        cost += distance(all_cities[tour[i]], all_cities[tour[i+1]]);
    }
    return cost;
}
