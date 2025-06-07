#include <iostream>
#include <vector>
#include <numeric>
#include <limits>
#include <algorithm>

double calculate_tour_cost(const std::vector<int>& tour, const std::vector<std::vector<double>>& dist_matrix) {
    double cost = 0.0;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
        cost += dist_matrix[tour[i]][tour[i+1]];
    }
    return cost;
}

/**
 * Solves the TSP by looking ahead one step with the greedy algorithm.
 */
std::vector<int> solve_optimized_lookahead_greedy(const std::vector<std::vector<double>>& dist_matrix) {
    const int num_cities = dist_matrix.size();
    if (num_cities == 0) {
        return {};
    }

    // For each city, create a sorted list of its neighbors by distance.
    std::vector<std::vector<int>> nearest_neighbors(num_cities);
    for (int i = 0; i < num_cities; ++i) {
        std::vector<std::pair<double, int>> neighbors;
        for (int j = 0; j < num_cities; ++j) {
            if (i != j) {
                neighbors.push_back({dist_matrix[i][j], j});
            }
        }
        std::sort(neighbors.begin(), neighbors.end());
        
        for (const auto& pair : neighbors) {
            nearest_neighbors[i].push_back(pair.second);
        }
    }

    std::vector<int> tour;
    tour.reserve(num_cities + 1);

    std::vector<bool> visited(num_cities, false);

    int current_city = 0;
    tour.push_back(current_city);
    visited[current_city] = true;

    while (tour.size() < num_cities) {
        int best_next_city = -1;
        double min_score = std::numeric_limits<double>::max();

        for (int candidate_city = 0; candidate_city < num_cities; ++candidate_city) {
            if (!visited[candidate_city]) {
                double primary_cost = dist_matrix[current_city][candidate_city];

                double lookahead_cost = 0.0;
                bool found_neighbor = false;
                for (int neighbor : nearest_neighbors[candidate_city]) {
                    if (!visited[neighbor]) {
                        lookahead_cost = dist_matrix[candidate_city][neighbor];
                        found_neighbor = true;
                        break; // Found the closest unvisited one, no need to look further
                    }
                }

                double total_score = primary_cost + lookahead_cost;

                //Greedy Choice
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

    tour.push_back(0);

    return tour;
}

int main() {
    std::vector<std::vector<double>> dist_matrix = {
        { 0, 10, 15, 20 },
        { 10, 0, 35, 25 },
        { 15, 35, 0, 30 },
        { 20, 25, 30, 0 }
    };

    std::cout << "Solving TSP" << std::endl;

    std::vector<int> tour = solve_optimized_lookahead_greedy(dist_matrix);
    double total_cost = calculate_tour_cost(tour, dist_matrix);

    std::cout << "Calculated Tour: ";
    for (size_t i = 0; i < tour.size(); ++i) {
        std::cout << tour[i] << (i == tour.size() - 1 ? "" : " -> ");
    }
    std::cout << std::endl;
    std::cout << "Total Tour Cost: " << total_cost << std::endl;

    return 0;
}
