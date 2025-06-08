#include "mst.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <climits>
#include <cmath>
#include <limits>

MSTTSPSolver::MSTTSPSolver() : m_n(0) {}

bool MSTTSPSolver::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "cerror: Could'nt open file " << filename << std::endl;
        return false;
    }
    if (!parseTSPLIBFormat(file)) {
        file.close();
        return false;
    }
    file.close();
    buildDistanceMatrix();
    return true;
}

std::pair<std::vector<int>, double> MSTTSPSolver::solve() {
    if (m_n <= 0) return {{}, 0.0};
    if (m_n == 1) return {{0, 0}, 0.0};

    std::vector<std::pair<int, int>> mst_edges = computeMST();
    buildMSTAdjacencyList(mst_edges);
    std::vector<bool> visited(m_n, false);
    std::vector<int> preorder;
    dfsPreorder(0, visited, preorder);
    std::vector<int> tour = shortcutTour(preorder);

    double total_cost = 0.0;
    for (size_t i = 0; i < tour.size() - 1; i++) {
        total_cost += m_dist[tour[i]][tour[i+1]];
    }
    return {tour, total_cost};
}


bool MSTTSPSolver::parseTSPLIBFormat(std::ifstream& file) {
    std::string line;
    bool reading_coords = false;
    while (std::getline(file, line)) {
        std::string key, value;
        std::istringstream iss(line);
        iss >> key;
        
        if (key == "NAME:") iss >> m_instance_name;
        else if (key == "DIMENSION:") {
            iss >> m_n;
            m_dist.resize(m_n, std::vector<double>(m_n, 0.0));
            m_coords.resize(m_n);
            m_mst_adj.resize(m_n);
        } else if (key == "NODE_COORD_SECTION") reading_coords = true;
        else if (key == "EOF") break;
        else if (reading_coords) {
            int id;
            double x, y;
            std::istringstream line_stream(line);
            if (line_stream >> id >> x >> y && id >= 1 && id <= m_n) {
                m_coords[id - 1] = {x, y};
            }
        }
    }
    return true;
}

double MSTTSPSolver::euclideanDistance(const std::pair<double, double>& p1, const std::pair<double, double>& p2) const {
    return sqrt(pow(p1.first - p2.first, 2) + pow(p1.second - p2.second, 2));
}

void MSTTSPSolver::buildDistanceMatrix() {
    for (int i = 0; i < m_n; i++) {
        for (int j = 0; j < m_n; j++) {
            m_dist[i][j] = (i == j) ? 0.0 : round(euclideanDistance(m_coords[i], m_coords[j]));
        }
    }
}

std::vector<std::pair<int, int>> MSTTSPSolver::computeMST() {
    std::vector<bool> in_mst(m_n, false);
    std::vector<double> key(m_n, std::numeric_limits<double>::max());
    std::vector<int> parent(m_n, -1);
    std::vector<std::pair<int, int>> mst_edges;
    
    key[0] = 0.0;

    for (int count = 0; count < m_n; count++) {
        int u = -1;
        for (int v = 0; v < m_n; v++) {
            if (!in_mst[v] && (u == -1 || key[v] < key[u])) {
                u = v;
            }
        }
        if(u == -1) break; // All remaining vertices are inaccessible
        in_mst[u] = true;
        if (parent[u] != -1) {
            mst_edges.push_back({parent[u], u});
        }
        for (int v = 0; v < m_n; v++) {
            if (!in_mst[v] && m_dist[u][v] < key[v]) {
                key[v] = m_dist[u][v];
                parent[v] = u;
            }
        }
    }
    return mst_edges;
}

void MSTTSPSolver::buildMSTAdjacencyList(const std::vector<std::pair<int, int>>& mst_edges) {
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
