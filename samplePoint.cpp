#include "samplePoint.h"

void samplePoint(vector<Point2d>& src_points, vector<Point2d>& dst_points, const int min_num, const double thresh_dist, vector<Point2d>& src_res, vector<Point2d>& dst_res) {
    default_random_engine e;
    vector<Point2d> src_out, dst_out;
    int count = src_points.size();
    while (true) {
        if (count <= 0) break;
        src_res.clear();
        dst_res.clear();
        src_out.clear();
        dst_out.clear();
        uniform_int_distribution<unsigned> u(0, src_points.size() - 1);
        int r = u(e);
        Point2d r_point = src_points[r];
        for (int i = 0; i < src_points.size(); ++i) {
            if (i == r) continue;
            double dist = sqrt(pow(r_point.x - src_points[i].x, 2) + pow(r_point.y - src_points[i].y, 2));
            if (dist <= thresh_dist) {
                src_res.push_back(src_points[i]);
                dst_res.push_back(dst_points[i]);
            }
            else {
                src_out.push_back(src_points[i]);
                dst_out.push_back(dst_points[i]);
            }
        }

        if (src_res.size() >= min_num) break;
        count--;
    }
    src_points.swap(src_out);
    dst_points.swap(dst_out);
}