#ifndef AANAP_SAMPLEPOINT_H
#define AANAP_SAMPLEPOINT_H

#include <opencv2/opencv.hpp>
#include <random>

using namespace std;
using namespace cv;

void samplePoint(vector<Point2d>& src_points, vector<Point2d>& dst_points, const int min_num, const double thresh_dist, vector<Point2d>& src_res, vector<Point2d>& dst_res);

#endif //AANAP_SAMPLEPOINT_H

