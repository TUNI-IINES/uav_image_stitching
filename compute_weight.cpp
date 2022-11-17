#include "compute_weight.h"

vector<double> compute_weight(const Point2d& p, const Point2d& K1, const Point2d& K2) {
    vector<double> r(2);
    if (p.x < K1.x) {
        r[0] = 0;
        r[1] = 1;
        return r;
    }

    if (p.x > K2.x) {
        r[0] = 1;
        r[1] = 0;
        return r;
    }

    double a = (p.x - K1.x) * (K2.x - K1.x);
    double b = (p.y - K1.y) * (K2.y - K1.y);
    double c = pow(K2.x - K1.x, 2) + pow(K2.y - K1.y, 2);

    r[0] = abs(a + b) / c;
    r[1] = 1 - r[0];
    return r;
}