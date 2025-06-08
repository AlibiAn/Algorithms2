#ifndef MST_H
#define MST_H

#include <vector>
#include <utility>
#include <string>
#include <fstream>

/*
    MST-based 2-approximation algorithm for TSP (CLRS Section 35.2.1)
    1) Computes a minimum spanning tree
    2) Takes a preorder tree walk of MST  
    3) Takes shortcutting to make a Hamiltonian cycle
    
    Approximation ratio: 2
    Time Complexity: O(V^2) where V is the number of vertices
*/

class MSTTSPSolver {
private: 
    // Core data members
    int m_n;                                          // # of vertices
    std::vector<std::vector<double>> m_dist;          // distance matrix  
    std::vector<std::pair<double, double>> m_coords;  // coordinates of vertices
    
    // Additional members for dataset support
    std::vector<std::vector<int>> m_mst_adj;          // MST adjacency list
    std::string m_instance_name;                      // Instance name from file
    std::string m_edge_weight_type;                   // Edge weight type (EUC_2D, etc.)
    
    // Core algorithm methods
    std::vector<std::pair<int, int>> computeMST();   // Prim's algorithm
    void buildMSTAdjacencyList(const std::vector<std::pair<int, int>>& mst_edges);
    void dfsPreorder(int vertex, std::vector<bool>& visited, std::vector<int>& preorder);
    std::vector<int> shortcutTour(const std::vector<int>& preorder);
    
    // Distance calculation methods
    double euclideanDistance(const std::pair<double, double>& p1, 
                           const std::pair<double, double>& p2) const;
    double euclideanDistanceRounded(const std::pair<double, double>& p1, 
                                  const std::pair<double, double>& p2) const;
    void buildDistanceMatrix();
    double getDistanceOnDemand(int i, int j) const;
    
    // File parsing methods
    bool parseTSPLIBFormat(std::ifstream& file);
    
public:
    explicit MSTTSPSolver(int vertices);
    MSTTSPSolver();
    
    // Basic functionality
    void setDistance(int i, int j, double distance);
    double getDistance(int i, int j) const;
    bool loadFromFile(const std::string& filename);
    std::pair<std::vector<int>, double> solve();
    double calculateTourCost(const std::vector<int>& tour);
    int getVertexCount() const;
    std::string getInstanceName() const;
};

#endif
