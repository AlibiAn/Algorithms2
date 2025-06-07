// Version 1 : Hybrid Look-Ahead Greedy + 2-Opt TSP Solver
//Main-ideas:
// 1.  Look-Ahead Greedy Constructor: Builds an initial tour by choosing the next city based on a
//     "look-ahead" score (cost to candidate + cost from candidate to its nearest neighbor).
// 2.  Multi-Start Randomization: Runs the constructor multiple times from different random starting
//     cities to avoid getting stuck in a poor solution due to a bad starting choice.
// 3.  2-Opt Local Search: After construction, a 2-Opt refinement phase is applied to iteratively
//     improve the tour by removing edge crossings.
//
//  Spotted after usage:
//    Onlt really works for small problems where a full N x N distance matrix can be held in memory.
//    Does not scale to large TSP instances from files(Kazakhstan dataset took too long to converge).
//    The constructor is O(N^2 log N) due to pre-computation, which is too slow for very large N.

#include <iostream>
#include <vector>
#include <numeric>
#include <limits>
#include <algorithm>
#include <random>

double calculate_tour_cost(const std::vector<int>& tour, const std::vector<std::vector<double>>& dist_matrix) {
    double cost = 0.0;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
        cost += dist_matrix[tour[i]][tour[i+1]];
    }
    return cost;
}

std::vector<int> solve_lookahead_greedy_constructor(
    const std::vector<std::vector<double>>& dist_matrix,
    const std::vector<std::vector<int>>& nearest_neighbors,
    int start_city
) {
    const int num_cities = dist_matrix.size();
    std::vector<int> tour;
    tour.reserve(num_cities + 1);
    std::vector<bool> visited(num_cities, false);

    int current_city = start_city;
    tour.push_back(current_city);
    visited[current_city] = true;

    while (tour.size() < num_cities) {
        int best_next_city = -1;
        double min_score = std::numeric_limits<double>::max();

        for (int candidate_city = 0; candidate_city < num_cities; ++candidate_city) {
            if (!visited[candidate_city]) {
                double primary_cost = dist_matrix[current_city][candidate_city];
                double lookahead_cost = 0.0;
                for (int neighbor : nearest_neighbors[candidate_city]) {
                    if (!visited[neighbor]) {
                        lookahead_cost = dist_matrix[candidate_city][neighbor];
                        break;
                    }
                }
                double total_score = primary_cost + lookahead_cost;
                if (total_score < min_score) {
                    min_score = total_score;
                    best_next_city = candidate_city;
                }
            }
        }
        if (best_next_city != -1) {
            current_city = best_next_city;
            tour.push_back(current_city);
            visited[current_city] = true;
        } else {
            break;
        }
    }
    tour.push_back(start_city);
    return tour;
}

void apply_two_opt(std::vector<int>& tour, const std::vector<std::vector<double>>& dist_matrix) {
    bool improved = true;
    int num_cities_in_tour = tour.size() - 1;
    while (improved) {
        improved = false;
        for (int i = 0; i < num_cities_in_tour - 1; ++i) {
            for (int j = i + 1; j < num_cities_in_tour; ++j) {
                double current_dist = dist_matrix[tour[i]][tour[i+1]] + dist_matrix[tour[j]][tour[j+1]];
                double new_dist = dist_matrix[tour[i]][tour[j]] + dist_matrix[tour[i+1]][tour[j+1]];
                if (new_dist < current_dist) {
                    std::reverse(tour.begin() + i + 1, tour.begin() + j + 1);
                    improved = true;
                }
            }
        }
    }
}

std::vector<int> solve_tsp_hybrid(const std::vector<std::vector<double>>& dist_matrix, int num_starts) {
    const int num_cities = dist_matrix.size();
    if (num_cities == 0) return {};
    if (num_starts <= 0) num_starts = 1;

    std::vector<std::vector<int>> nearest_neighbors(num_cities);
    for (int i = 0; i < num_cities; ++i) {
        std::vector<std::pair<double, int>> neighbors;
        for (int j = 0; j < num_cities; ++j) {
            if (i != j) neighbors.push_back({dist_matrix[i][j], j});
        }
        std::sort(neighbors.begin(), neighbors.end());
        for (const auto& pair : neighbors) {
            nearest_neighbors[i].push_back(pair.second);
        }
    }

    std::vector<int> best_tour;
    double min_cost = std::numeric_limits<double>::max();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0, num_cities - 1);
    
    for (int i = 0; i < num_starts; ++i) {
        int start_city = distrib(gen);
        std::vector<int> current_tour = solve_lookahead_greedy_constructor(dist_matrix, nearest_neighbors, start_city);
        apply_two_opt(current_tour, dist_matrix);
        double current_cost = calculate_tour_cost(current_tour, dist_matrix);
        if (current_cost < min_cost) {
            min_cost = current_cost;
            best_tour = current_tour;
        }
    }
    return best_tour;
}

int main() {
    std::vector<std::vector<double>> dist_matrix = {
        { 0, 10, 15, 20 },
        { 10, 0, 35, 25 },
        { 15, 35, 0, 30 },
        { 20, 25, 30, 0 }
    };
    int num_starts = dist_matrix.size();

    std::cout << "Version 1: " << std::endl;
    std::vector<int> tour = solve_tsp_hybrid(dist_matrix, num_starts);
    double total_cost = calculate_tour_cost(tour, dist_matrix);

    std::cout << "Best tour found: ";
    for (size_t i = 0; i < tour.size(); ++i) {
        std::cout << tour[i] << (i == tour.size() - 1 ? "" : " -> ");
    }
    std::cout << std::endl;
    std::cout << "Total Tour Cost: " << total_cost << std::endl;

    return 0;
}
