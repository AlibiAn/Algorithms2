#ifndef HELD_KARP_H
#define HELD_KARP_H

#include <vector>
#include <utility>
#include <string>
#include <fstream>

class HeldKarpTSPSolver {
private:
    int m_n;
    std::vector<std::vector<double>> m_dist;
    std::vector<std::pair<double, double>> m_coords;
    std::string m_instance_name;

    std::vector<std::vector<double>> m_dp;
    std::vector<std::vector<int>> m_parent;

    double euclideanDistance(const std::pair<double, double>& p1, const std::pair<double, double>& p2) const;
    void buildDistanceMatrix();
    bool parseTSPLIBFormat(std::ifstream& file);
    
    void solveDP();
    std::vector<int> reconstructPath(int last_vertex);

public:
    HeldKarpTSPSolver();
    
    bool loadFromFile(const std::string& filename);
    std::pair<std::vector<int>, double> solve();
    int getVertexCount() const;
    std::string getInstanceName() const;
};

#endif
