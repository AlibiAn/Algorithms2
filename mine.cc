
// Version 3: Clustering with Closest Node Pair Stitching
//
//More Changes from V2:
// 1.  Upgraded Tour Stitching Logic: This is the key improvement. The Combining phase was
//     completely reworked to be better.
// 2.  New dedicated function was created to find
//     the absolute shortest "bridge" between two sub-tours by checking every possible pair of nodes.
// 3.  Merging: Instead of connecting a sub-tour to the abstract centroid of its neighbor,
//     the algorithm now connects it to the actual nearest node.
//
// New Insights:
//  Slightly better than V2 at non centroid data like Kazakhstan map
// 

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
std::pair<int, int> find_closest_connection(const std::vector<int>& tour1, const std::vector<int>& tour2, const std::vector<City>& all_cities);
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

    std::cout << "Version 3: " << std::endl;
    std::cout << "Successfully parsed " << cities.size() << " cities." << std::endl;
    std::cout << "Running with K=" << K << " and " << sub_problem_starts << " starts for sub-problems." << std::endl;

    std::vector<int> final_tour = solve_tsp_cluster_hybrid(cities, K, sub_problem_starts);
    double final_cost = calculate_tour_cost(final_tour, cities);

    std::cout << "\n" << " Best Result Found\n" << std::endl;
    std::cout << "Final Tour Cost: " << final_cost << std::endl;

    return 0;
}

// --- Main Orchestrator with NEW Stitching Logic ---
std::vector<int> solve_tsp_cluster_hybrid(const std::vector<City>& cities, int K, int sub_problem_starts) {
    std::cout << "\nClustering cities" << std::endl;
    auto [clusters_of_global_ids, centroids] = kmeans_plus_plus(cities, K, 100);

    std::cout << "\nSolving TSP for each cluster" << std::endl;
    std::vector<std::vector<int>> sub_tours_global(K);
    for (int i = 0; i < K; ++i) {
        if (clusters_of_global_ids[i].empty()) continue;
        std::vector<City> sub_problem_cities;
        for (int global_id : clusters_of_global_ids[i]) {
            sub_problem_cities.push_back(cities[global_id]);
        }
        sub_tours_global[i] = solve_sub_problem(sub_problem_cities, sub_problem_starts, cities);
    }
    
    std::cout << "\nFinding optimal cluster sequence" << std::endl;
    std::vector<int> centroid_tour_local = solve_sub_problem(centroids, K * 2, centroids);
    if (centroid_tour_local.empty()) return {};
    centroid_tour_local.pop_back();

    std::cout << "\nStitching sub-tours with the Closest Node pair" << std::endl;
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

    return final_tour;
}

std::pair<int, int> find_closest_connection(const std::vector<int>& tour1, const std::vector<int>& tour2, const std::vector<City>& all_cities) {
    int best_node1 = -1, best_node2 = -1;
    double min_dist = std::numeric_limits<double>::max();

    if (tour1.empty() || tour2.empty()) return {-1, -1};

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
    return {best_node1, best_node2};
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
