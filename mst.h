#ifndef MST_H
#define MST_H

#include <vector>
#include <utility>
#include <string>
#include <fstream>

class MSTTSPSolver {
private: 
    int m_n;
    std::vector<std::pair<double, double>> m_coords;
    std::vector<std::vector<int>> m_mst_adj;
    std::string m_instance_name;

    std::vector<std::pair<int, int>> computeMST();
    void buildMSTAdjacencyList(const std::vector<std::pair<int, int>>& mst_edges);
    void dfsPreorder(int vertex, std::vector<bool>& visited, std::vector<int>& preorder);
    std::vector<int> shortcutTour(const std::vector<int>& preorder);
    
    double euclideanDistance(const std::pair<double, double>& p1, const std::pair<double, double>& p2) const;
    bool parseTSPLIBFormat(std::ifstream& file);
    
    double getDistance(int i, int j) const;

public:
    MSTTSPSolver();
    
    bool loadFromFile(const std::string& filename);
    std::pair<std::vector<int>, double> solve();
    int getVertexCount() const;
    std::string getInstanceName() const;
};

#endif
