#ifndef AANAP_COMPUTE_WEIGHT_H
#define AANAP_COMPUTE_WEIGHT_H

#include <iostream>
#include "opencv2/opencv.hpp"

using namespace std;
using namespace cv;

vector<double> compute_weight(const Point2d& p, const Point2d& K1, const Point2d& K2);

#endif //AANAP_COMPUTE_WEIGHT_H

