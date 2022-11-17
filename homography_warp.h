#ifndef AANAP_HOMOGRAPHY_WARP_H
#define AANAP_HOMOGRAPHY_WARP_H

#include <opencv2/opencv.hpp>

using namespace cv;
using namespace std;

Point2d homography_warp(const Mat& src, const Mat& H, const Point2d offset, const Size s, Mat& dst);

Mat homography_linearization(const Mat& H, const Point& center_point, const vector<Point2d>& anchor_points);

#endif //AANAP_HOMOGRAPHY_WARP_H
