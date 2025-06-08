#include "mst.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <climits>
#include <cmath>
#include <limits>
#include <queue>

MSTTSPSolver::MSTTSPSolver() : m_n(0) {}

bool MSTTSPSolver::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    if (!parseTSPLIBFormat(file)) {
        file.close();
        return false;
    }
    file.close();
    return true;
}

std::pair<std::vector<int>, double> MSTTSPSolver::solve() {
    if (m_n <= 0) return {{}, 0.0};
    if (m_n == 1) return {{0, 0}, 0.0};

    auto mst_edges = computeMST();
    buildMSTAdjacencyList(mst_edges);
    std::vector<bool> visited(m_n, false);
    std::vector<int> preorder;
    dfsPreorder(0, visited, preorder);
    auto tour = shortcutTour(preorder);

    double total_cost = 0.0;
    for (size_t i = 0; i < tour.size() - 1; i++) {
        total_cost += getDistance(tour[i], tour[i+1]);
    }
    return {tour, total_cost};
}

bool MSTTSPSolver::parseTSPLIBFormat(std::ifstream& file) {
    std::string line;
    bool reading_coords = false;
    while (std::getline(file, line)) {
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (line.empty()) continue;
        if (line.rfind("NAME", 0) == 0) {
            size_t colon_pos = line.find(':');
            if (colon_pos != std::string::npos) {
                m_instance_name = line.substr(colon_pos + 1);
                m_instance_name.erase(0, m_instance_name.find_first_not_of(" \t"));
            }
        } else if (line.rfind("DIMENSION", 0) == 0) {
            size_t colon_pos = line.find(':');
            std::string dim_str = (colon_pos != std::string::npos) ? line.substr(colon_pos + 1) : line.substr(9);
            dim_str.erase(0, dim_str.find_first_not_of(" \t"));
            try { m_n = std::stoi(dim_str); } catch (...) { return false; }
            if (m_n > 0) {
                m_coords.resize(m_n);
                m_mst_adj.resize(m_n);
            }
        } else if (line == "NODE_COORD_SECTION") {
            reading_coords = true;
        } else if (line == "EOF") {
            break;
        } else if (reading_coords) {
            std::istringstream line_stream(line);
            int id;
            double x, y;
            if (line_stream >> id >> x >> y && id >= 1 && id <= m_n) {
                m_coords[id - 1] = {x, y};
            }
        }
    }
    return m_n > 0;
}

double MSTTSPSolver::euclideanDistance(const std::pair<double, double>& p1, const std::pair<double, double>& p2) const {
    return sqrt(pow(p1.first - p2.first, 2) + pow(p1.second - p2.second, 2));
}

double MSTTSPSolver::getDistance(int i, int j) const {
    if (i == j) return 0.0;
    return round(euclideanDistance(m_coords[i], m_coords[j]));
}

std::vector<std::pair<int, int>> MSTTSPSolver::computeMST() {
    std::vector<double> key(m_n, std::numeric_limits<double>::max());
    std::vector<int> parent(m_n, -1);
    std::vector<bool> in_mst(m_n, false);
    std::vector<std::pair<int, int>> mst_edges;
    
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> pq;
    
    key[0] = 0.0;
    pq.push({0.0, 0});
    
    while (!pq.empty()) {
        int u = pq.top().second;
        pq.pop();
        if (in_mst[u]) continue;
        
        in_mst[u] = true;
        if (parent[u] != -1) {
            mst_edges.push_back({parent[u], u});
        }
        
        for (int v = 0; v < m_n; v++) {
            if (!in_mst[v]) {
                double weight = getDistance(u, v);
                if (weight < key[v]) {
                    key[v] = weight;
                    parent[v] = u;
                    pq.push({weight, v});
                }
            }
        }
    }
    return mst_edges;
}

void MSTTSPSolver::buildMSTAdjacencyList(const std::vector<std::pair<int, int>>& mst_edges) {
    for (int i = 0; i < m_n; i++) m_mst_adj[i].clear();
    for (const auto& edge : mst_edges) {
        m_mst_adj[edge.first].push_back(edge.second);
        m_mst_adj[edge.second].push_back(edge.first);
    }
}

void MSTTSPSolver::dfsPreorder(int vertex, std::vector<bool>& visited, std::vector<int>& preorder) {
    visited[vertex] = true;
    preorder.push_back(vertex);
    for (int neighbor : m_mst_adj[vertex]) {
        if (!visited[neighbor]) {
            dfsPreorder(neighbor, visited, preorder);
        }
    }
}

std::vector<int> MSTTSPSolver::shortcutTour(const std::vector<int>& preorder) {
    std::vector<bool> visited(m_n, false);
    std::vector<int> tour;
    for (int vertex : preorder) {
        if (!visited[vertex]) {
            visited[vertex] = true;
            tour.push_back(vertex);
        }
    }
    if (!tour.empty()) {
        tour.push_back(tour.front());
    }
    return tour;
}

int MSTTSPSolver::getVertexCount() const { return m_n; }
std::string MSTTSPSolver::getInstanceName() const { return m_instance_name; }
