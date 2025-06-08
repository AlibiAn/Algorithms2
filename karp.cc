#include "karp.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <climits>
#include <cmath>
#include <limits>

HeldKarpTSPSolver::HeldKarpTSPSolver() : m_n(0) {}

bool HeldKarpTSPSolver::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return false;
    }
    if (!parseTSPLIBFormat(file)) {
        file.close();
        return false;
    }
    file.close();

    if (m_n > 0) {
        if (m_n > 25) {
            std::cout << "Skipping held-karp: " << m_n << " vertices is too large(crashes Ubuntu)" << std::endl;
            return true;
        }
        int num_subsets = 1 << m_n;
        m_dp.assign(num_subsets, std::vector<double>(m_n, std::numeric_limits<double>::max()));
        m_parent.assign(num_subsets, std::vector<int>(m_n, -1));
        buildDistanceMatrix();
    }
    return true;
}

std::pair<std::vector<int>, double> HeldKarpTSPSolver::solve() {
    if (m_n <= 0) return {{}, -1.0};
    if (m_n > 25) return {{}, -1.0};
    if (m_n == 1) return {{0, 0}, 0.0};

    solveDP();

    int final_mask = (1 << m_n) - 1;
    double min_cost = std::numeric_limits<double>::max();
    int last_vertex = -1;

    for (int i = 1; i < m_n; ++i) {
        if (m_dp[final_mask][i] != std::numeric_limits<double>::max()) {
            double current_cost = m_dp[final_mask][i] + m_dist[i][0];
            if (current_cost < min_cost) {
                min_cost = current_cost;
                last_vertex = i;
            }
        }
    }

    if (last_vertex == -1) return {{}, -1.0}; 
    std::vector<int> tour = reconstructPath(last_vertex);
    return {tour, min_cost};
}

void HeldKarpTSPSolver::solveDP() {
    m_dp[1][0] = 0.0; 

    for (int mask = 1; mask < (1 << m_n); mask += 2) {
        for (int u = 0; u < m_n; ++u) {
            if (!(mask & (1 << u))) continue; 
            if (m_dp[mask][u] == std::numeric_limits<double>::max()) continue;

            for (int v = 0; v < m_n; ++v) {
                if (!(mask & (1 << v))) {
                    int next_mask = mask | (1 << v);
                    double new_cost = m_dp[mask][u] + m_dist[u][v];
                    if (new_cost < m_dp[next_mask][v]) {
                        m_dp[next_mask][v] = new_cost;
                        m_parent[next_mask][v] = u;
                    }
                }
            }
        }
    }
}

std::vector<int> HeldKarpTSPSolver::reconstructPath(int last_vertex) {
    std::vector<int> path;
    int current_mask = (1 << m_n) - 1;
    int current_vertex = last_vertex;

    while (current_vertex != -1) {
        path.push_back(current_vertex);
        int prev_vertex = m_parent[current_mask][current_vertex];
        current_mask ^= (1 << current_vertex);
        current_vertex = prev_vertex;
    }
    std::reverse(path.begin(), path.end());
    path.push_back(0);
    return path;
}

bool HeldKarpTSPSolver::parseTSPLIBFormat(std::ifstream& file) {
    std::string line;
    bool reading_coords = false;
    while (std::getline(file, line)) {
        std::string key;
        std::istringstream iss(line);
        iss >> key;
        if (key == "NAME:") iss >> m_instance_name;
        else if (key == "DIMENSION:") iss >> m_n;
        else if (key == "NODE_COORD_SECTION") reading_coords = true;
        else if (key == "EOF") break;
        else if (reading_coords) {
            int id; double x, y;
            std::istringstream line_stream(line);
            if (line_stream >> id >> x >> y && id >= 1 && id <= m_n) {
                if(m_coords.empty()) m_coords.resize(m_n);
                m_coords[id - 1] = {x, y};
            }
        }
    }
    return m_n > 0;
}

double HeldKarpTSPSolver::euclideanDistance(const std::pair<double, double>& p1, const std::pair<double, double>& p2) const {
    return sqrt(pow(p1.first - p2.first, 2) + pow(p1.second - p2.second, 2));
}

void HeldKarpTSPSolver::buildDistanceMatrix() {
    for (int i = 0; i < m_n; ++i) {
        for (int j = 0; j < m_n; ++j) {
            m_dist[i][j] = (i == j) ? 0.0 : round(euclideanDistance(m_coords[i], m_coords[j]));
        }
    }
}

int HeldKarpTSPSolver::getVertexCount() const { return m_n; }
std::string HeldKarpTSPSolver::getInstanceName() const { return m_instance_name; }
