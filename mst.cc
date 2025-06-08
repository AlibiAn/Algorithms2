#include "mst.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <climits>
#include <cmath>
#include <cctype>
#include <limits>

MSTTSPSolver::MSTTSPSolver(int vertices) : m_n(vertices), m_edge_weight_type("EUC_2D") {
    if (m_n > 0) {
        // For large instances, avoid precomputing full distance matrix
        if (m_n <= 10000) {
            m_dist.resize(m_n, std::vector<double>(m_n, 0.0));
        }
        m_mst_adj.resize(m_n);
        m_coords.resize(m_n);
    }
}

MSTTSPSolver::MSTTSPSolver() : m_n(0), m_edge_weight_type("EUC_2D") {}

void MSTTSPSolver::setDistance(int i, int j, double distance) {
    if (i >= 0 && i < m_n && j >= 0 && j < m_n && !m_dist.empty()) {
        m_dist[i][j] = m_dist[j][i] = distance;
    }
}

double MSTTSPSolver::euclideanDistance(const std::pair<double, double>& p1, 
                                     const std::pair<double, double>& p2) const {
    double dx = p1.first - p2.first;
    double dy = p1.second - p2.second;
    return sqrt(dx * dx + dy * dy);
}

double MSTTSPSolver::euclideanDistanceRounded(const std::pair<double, double>& p1, 
                                            const std::pair<double, double>& p2) const {
    double distance = euclideanDistance(p1, p2);
    return round(distance);
}

double MSTTSPSolver::getDistanceOnDemand(int i, int j) const {
    if (i == j) return 0.0;
    
    // If we have precomputed distance matrix, use it
    if (!m_dist.empty() && i < (int)m_dist.size() && j < (int)m_dist[i].size()) {
        return m_dist[i][j];
    }
    
    // Otherwise compute on demand from coordinates
    if (i < (int)m_coords.size() && j < (int)m_coords.size()) {
        if (m_edge_weight_type == "EUC_2D") {
            return euclideanDistanceRounded(m_coords[i], m_coords[j]);
        } else {
            return euclideanDistance(m_coords[i], m_coords[j]);
        }
    }
    
    return std::numeric_limits<double>::max();
}

void MSTTSPSolver::buildDistanceMatrix() {
    // Only build full matrix for smaller instances
    if (m_n <= 10000) {
        m_dist.resize(m_n, std::vector<double>(m_n, 0.0));
        for (int i = 0; i < m_n; i++) {
            for (int j = 0; j < m_n; j++) {
                if (i == j) {
                    m_dist[i][j] = 0.0;
                } else {
                    if (m_edge_weight_type == "EUC_2D") {
                        m_dist[i][j] = euclideanDistanceRounded(m_coords[i], m_coords[j]);
                    } else {
                        m_dist[i][j] = euclideanDistance(m_coords[i], m_coords[j]);
                    }
                }
            }
        }
    }
}


bool MSTTSPSolver::parseTSPLIBFormat(std::ifstream& file) {
    std::string line;
    bool reading_coords = false;
    
    while (std::getline(file, line)) {
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        
        if (line.empty()) continue;
        
        if (line.find("NAME") == 0) {
            size_t colon_pos = line.find(':');
            if (colon_pos != std::string::npos) {
                m_instance_name = line.substr(colon_pos + 1);
                m_instance_name.erase(0, m_instance_name.find_first_not_of(" \t"));
                m_instance_name.erase(m_instance_name.find_last_not_of(" \t") + 1);
            }
            
        } else if (line.find("DIMENSION") == 0) {
            std::istringstream iss(line);
            std::string token;
            iss >> token;
            if (line.find(":") != std::string::npos) {
                iss >> token;
            }
            iss >> m_n;
            
            std::cout << "Loading instance: " << m_instance_name 
                      << " with " << m_n << " vertices" << std::endl;
            
            // Resize data structures
            if (m_n <= 10000) {
                m_dist.resize(m_n, std::vector<double>(m_n, 0.0));
            }
            m_mst_adj.resize(m_n);
            m_coords.resize(m_n);
            
        } else if (line.find("EDGE_WEIGHT_TYPE") == 0) {
            std::istringstream iss(line);
            std::string token;
            iss >> token;
            if (line.find(":") != std::string::npos) {
                iss >> token;
            }
            iss >> m_edge_weight_type;
            
        } else if (line == "NODE_COORD_SECTION") {
            reading_coords = true;            
        } else if (line == "EOF") {
            break;
            
        } else if (reading_coords && !line.empty()) {
            std::istringstream iss(line);
            int id;
            double x, y;
            
            if (iss >> id >> x >> y) {
                if (id >= 1 && id <= m_n) {
                    m_coords[id - 1] = {x, y};
                }
            }
        }
    }
    
    return true;
}

bool MSTTSPSolver::loadFromFile(const std::string& filename) {
    std::cout << "Loading file: " << filename << std::endl;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Cannot open file " << filename << std::endl;
        return false;
    }
    
    bool success = parseTSPLIBFormat(file);
    file.close();
    
    if (!success) {
        std::cerr << "Error parsing the file" << std::endl;
        return false;
    }
    
    std::cout << "Building distance matrix..." << std::endl;
    buildDistanceMatrix();
    
    std::cout << "Successfully loaded " << m_instance_name 
              << " (" << m_n << " vertices)" << std::endl;
    return true;
}

std::vector<std::pair<int, int>> MSTTSPSolver::computeMST() {
    std::vector<bool> in_mst(m_n, false);
    std::vector<double> key(m_n, std::numeric_limits<double>::max());
    std::vector<int> parent(m_n, -1);
    std::vector<std::pair<int, int>> mst_edges;
    
    key[0] = 0.0;
    
    int progress_step = m_n / 10;
    if (progress_step == 0) progress_step = 1;
    
    for (int count = 0; count < m_n; count++) {
        if (m_n > 1000 && count % progress_step == 0) {
            std::cout << "MST Progress: " << (count * 100 / m_n) << "%" << std::endl;
        }
        
        int u = -1;
        for (int v = 0; v < m_n; v++) {
            if (!in_mst[v] && (u == -1 || key[v] < key[u])) {
                u = v;
            }
        }
        
        in_mst[u] = true;
        
        if (parent[u] != -1) {
            mst_edges.push_back({parent[u], u});
        }
        
        // Update key values of adjacent vertices
        for (int v = 0; v < m_n; v++) {
            if (!in_mst[v]) {
                double weight = getDistanceOnDemand(u, v);
                if (weight < key[v]) {
                    key[v] = weight;
                    parent[v] = u;
                }
            }
        }
    }
    
    if (m_n > 1000) {
        std::cout << "MST computation completed." << std::endl;
    }
    
    return mst_edges;
}

void MSTTSPSolver::buildMSTAdjacencyList(const std::vector<std::pair<int, int>>& mst_edges) {
    // Clear existing adjacency list
    for (int i = 0; i < m_n; i++) {
        m_mst_adj[i].clear();
    }
    
    for (const auto& edge : mst_edges) {
        m_mst_adj[edge.first].push_back(edge.second);
        m_mst_adj[edge.second].push_back(edge.first);
    }
}

void MSTTSPSolver::dfsPreorder(int vertex, std::vector<bool>& visited, std::vector<int>& preorder) {
    visited[vertex] = true;
    preorder.push_back(vertex);
    
    // Visit all unvisited neighbors
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
        tour.push_back(tour[0]);
    }
    
    return tour;
}

std::pair<std::vector<int>, double> MSTTSPSolver::solve() {
    if (m_n <= 0) {
        return {{}, 0.0};
    }
    
    if (m_n == 1) {
        return {{0, 0}, 0.0};
    }
    
    std::cout << "\nSolving TSP using MST" << std::endl;
    
    std::cout << "Computing mst" << std::endl;
    std::vector<std::pair<int, int>> mst_edges = computeMST();
    
    std::cout << "Making the adjacency list of MST" << std::endl;
    buildMSTAdjacencyList(mst_edges);
    
    std::cout << "DFS (preorder) traversal" << std::endl;
    std::vector<bool> visited(m_n, false);
    std::vector<int> preorder;
    dfsPreorder(0, visited, preorder);
    
    std::cout << "Apply shortcuts" << std::endl;
    std::vector<int> tour = shortcutTour(preorder);
    
    std::cout << "Calculating tour cost" << std::endl;
    double total_cost = calculateTourCost(tour);
    
    return {tour, total_cost};
}

double MSTTSPSolver::calculateTourCost(const std::vector<int>& tour) {
    if (tour.size() < 2) {
        return 0.0;
    }
    
    double total_cost = 0.0;
    for (size_t i = 0; i < tour.size() - 1; i++) {
        total_cost += getDistanceOnDemand(tour[i], tour[i + 1]);
    }
    
    return total_cost;
}

int MSTTSPSolver::getVertexCount() const {
    return m_n;
}

std::string MSTTSPSolver::getInstanceName() const {
    return m_instance_name;
}

double MSTTSPSolver::getDistance(int i, int j) const {
    return getDistanceOnDemand(i, j);
}
