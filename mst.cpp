#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

struct Coordinate {
    double x, y;
};

double distance(Coordinate a, Coordinate b) {
    double d_x = a.x - b.x;
    double d_y = a.y - b.y;
    return sqrt((d_x*d_x) + (d_y * d_y));
}