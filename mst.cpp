#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>

using namespace std;

struct Coordinate {
    double x, y;
};

vector<Coordinate> accessTSP(int& n, string filename) {
    ifstream file(filename);
    string line;

    while (getline(file, line) && line.find("DIMENSION") == string::npos);
    sscanf(line.c_str(), "DIMENSION : %d", &n);

    while (getline(file, line) && line.find("NODE_COORD_SECTION") == string::npos);
    vector<Coordinate> coordinates(n);
    
    for (int j = 0; j < n; j++) {
        int id;
        file >> id >> points[i].x >> points[i].y;
    }

    return coordinates;

}

double distance(Coordinate a, Coordinate b) {
    double d_x = a.x - b.x;
    double d_y = a.y - b.y;
    return sqrt((d_x*d_x) + (d_y * d_y));
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cout << "Used" << argv[0] << " <tsp_file>" << endl;
        return -1;
    }

    return 0;
}