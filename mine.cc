// Version 2 (V2): Clustering-Based Divide and Conquer with Centroid Stitching
//
//Core Changes from V1:
// 1.  Scalability for Large Datasets: The entire architecture was changed to handle large problems.
//      Immediate Distance Calculation, that avoids creating a huge N x N distance matrix.
// 2.  Divide and Conquer Strategy: Implements a basic three-phase approach.
//     -   Phase 1 (Divide): Uses K-Means++ to partition cities into K clusters.
//     -   Phase 2 (Conquer): A sub-problem solver finds tours within each cluster.
//     -   Phase 3 (Combine): Finds the optimal cluster sequence on the centroids.
//
// More Problems:
//   This stitching logic is fast but not perfect. It performs poorly on datasets with non-spherical
//     clusters (like Kazakhstan) because a centroid is a poor representative of an irregular shape.
//      But Mona_Lisa and the rest of the datasets performed at a good level.

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <numeric>
#include <limits>
#include <algorithm>
#include <random>
#include <map>

struct City { int id; double x, y; };

bool parse_tsplib_file(const std::string& filename, std::vector<City>& cities);
double distance(const City& c1, const City& c2);
double calculate_tour_cost(const std::vector<int>& tour, const std::vector<City>& all_cities);
std::vector<int> solve_tsp_cluster_hybrid(const std::vector<City>& cities, int K, int sub_problem_starts);
std::pair<std::vector<std::vector<int>>, std::vector<City>> kmeans_plus_plus(const std::vector<City>& cities, int K, int max_iter);
std::vector<int> solve_sub_problem(const std::vector<City>& sub_cities, int num_starts, const std::vector<City>& all_cities_ref);


int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <path_to_tsp_file> <num_clusters> [sub_problem_starts]" << std::endl;
        return 1;
    }
    std::string filename = argv[1];
    int K = std::stoi(argv[2]);
    int sub_problem_starts = (argc > 3) ? std::stoi(argv[3]) : 10;

    std::vector<City> cities;
    if (!parse_tsplib_file(filename, cities)) return 1;

    std::cout << "Version 2: " << std::endl;
    std::cout << "Successfully parsed " << cities.size() << " cities." << std::endl;
    std::cout << "Running with K=" << K << " and " << sub_problem_starts << " starts for sub-problems." << std::endl;

    std::vector<int> final_tour = solve_tsp_cluster_hybrid(cities, K, sub_problem_starts);
    double final_cost = calculate_tour_cost(final_tour, cities);

    std::cout << "\n" << "           Best Result Found\n"  << std::endl;
    std::cout << "Final Tour Cost: " << final_cost << std::endl;

    return 0;
}


std::vector<int> solve_tsp_cluster_hybrid(const std::vector<City>& cities, int K, int sub_problem_starts) {
    std::cout << "\n Clustering Cities" << std::endl;
    auto [clusters_of_global_ids, centroids] = kmeans_plus_plus(cities, K, 100);

    std::cout << "\nSolving TSP for each cluster" << std::endl;
    std::vector<std::vector<int>> sub_tours_global(K);
    for (int i = 0; i < K; ++i) {
        if (clusters_of_global_ids[i].empty()) continue;
        std::vector<City> sub_problem_cities;
        for (int global_id : clusters_of_global_ids[i]) sub_problem_cities.push_back(cities[global_id]);
        sub_tours_global[i] = solve_sub_problem(sub_problem_cities, sub_problem_starts, cities);
    }
    
    std::cout << "\nFinding optimal cluster sequence" << std::endl;
    std::vector<int> centroid_tour_local = solve_sub_problem(centroids, K * 2, centroids);
    if (centroid_tour_local.empty()) return {};
    centroid_tour_local.pop_back(); // Get sequence from tour

    std::cout << "\n Stitching sub-tours with Centroid logic" << std::endl;
    std::vector<int> final_tour;
    
    for (size_t i = 0; i < centroid_tour_local.size(); ++i) {
        int current_cluster_idx = centroid_tour_local[i];
        if (sub_tours_global[current_cluster_idx].empty()) continue;
        
        int prev_cluster_idx = (i == 0) ? centroid_tour_local.back() : centroid_tour_local[i - 1];
        const City& prev_centroid = centroids[prev_cluster_idx];
        std::vector<int> current_sub_tour = sub_tours_global[current_cluster_idx];
        current_sub_tour.pop_back();

        int best_entry_node_global_id = -1;
        double min_entry_dist = std::numeric_limits<double>::max();
        for (int global_id : current_sub_tour) {
            double d = distance(cities[global_id], prev_centroid);
            if (d < min_entry_dist) {
                min_entry_dist = d;
                best_entry_node_global_id = global_id;
            }
        }
        
        auto it = std::find(current_sub_tour.begin(), current_sub_tour.end(), best_entry_node_global_id);
        if (it != current_sub_tour.end()) {
            std::rotate(current_sub_tour.begin(), it, current_sub_tour.end());
        }
        
        final_tour.insert(final_tour.end(), current_sub_tour.begin(), current_sub_tour.end());
    }

    if (!final_tour.empty()) final_tour.push_back(final_tour.front());

    return final_tour;
}

void apply_two_opt_local(std::vector<int>& local_tour, const std::vector<City>& sub_cities) {
    bool improved = true;
    int n = local_tour.size() - 1;
    while (improved) {
        improved = false;
        for (int i = 0; i < n - 1; ++i) {
            for (int j = i + 1; j < n; ++j) {
                const City& c_i = sub_cities[local_tour[i]];
                const City& c_i1 = sub_cities[local_tour[i+1]];
                const City& c_j = sub_cities[local_tour[j]];
                const City& c_j1 = sub_cities[local_tour[j+1]];
                if (distance(c_i, c_j) + distance(c_i1, c_j1) < distance(c_i, c_i1) + distance(c_j, c_j1)) {
                    std::reverse(local_tour.begin() + i + 1, local_tour.begin() + j + 1);
                    improved = true;
                }
            }
        }
    }
}
std::vector<int> solve_sub_problem(const std::vector<City>& sub_cities, int num_starts, const std::vector<City>& all_cities_ref) {
    const int num_sub_cities = sub_cities.size();
    if (num_sub_cities == 0) return {};
    if (num_sub_cities == 1) return {sub_cities[0].id, sub_cities[0].id};

    std::vector<int> best_tour_local;
    double min_cost = std::numeric_limits<double>::max();
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0, num_sub_cities - 1);

    for (int i = 0; i < num_starts; ++i) {
        int start_node_local = distrib(gen);
        std::vector<int> current_tour_local;
        std::vector<bool> visited(num_sub_cities, false);
        current_tour_local.push_back(start_node_local);
        visited[start_node_local] = true;

        while(current_tour_local.size() < num_sub_cities){
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
            if(best_next_local != -1) { current_tour_local.push_back(best_next_local); visited[best_next_local] = true; } else break;
        }
        current_tour_local.push_back(start_node_local);
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
std::pair<std::vector<std::vector<int>>, std::vector<City>> kmeans_plus_plus(const std::vector<City>& cities, int K, int max_iter) {
    const int N = cities.size();
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
            for (int k = 0; k < centroids.size(); ++k) {
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
    if (!file.is_open()) { std::cerr << "Error: Could not open file " << filename << std::endl; return false; }
    std::string line;
    bool in_coord_section = false;
    while (std::getline(file, line)) {
        if (line.find("NODE_COORD_SECTION") != std::string::npos) { in_coord_section = true; continue; }
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
