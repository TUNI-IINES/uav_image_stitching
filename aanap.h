#ifndef AANAP_APAP_H
#define AANAP_APAP_H

#include <opencv2/stitching/detail/warpers.hpp>
#include <opencv2/stitching/detail/matchers.hpp>


using namespace cv;
using namespace cv::detail;
using namespace std;

typedef double apap_float;

class AANAPWarper
{
public:
    AANAPWarper();
    ~AANAPWarper();

    int buildMaps(cv::Mat src_img, cv::detail::ImageFeatures src_features, ImageFeatures dst_features,
        MatchesInfo matches_info, Mat& xmap, Mat& ymap, Point& corner, Mat dst_img);
    int buildMaps(vector<Mat> imgs, vector<ImageFeatures> features,
        vector<MatchesInfo> pairwise_matches,
        vector<Mat>& xmaps, vector<Mat>& ymaps, vector<Point>& corners);

protected:

private:
    int CellDLT(int offset_x, int offset_y, ImageFeatures src_features,
        ImageFeatures dst_features, MatchesInfo matches_info, Mat_<apap_float>& H, Mat_<apap_float>& H_r, const Mat& dst_img);
    void InitA(ImageFeatures src_features, ImageFeatures dst_features,
        MatchesInfo matches_info);

    void FindS(ImageFeatures src_features, ImageFeatures dst_features,
        MatchesInfo matches_info, const Mat src);

    void testh(ImageFeatures src_features, ImageFeatures dst_features,
        MatchesInfo matches_info, Mat h);

    void BuildCellMap(int x1, int y1, int x2, int y2, Mat H, Point corner, Mat& xmap, Mat& ymap);

    inline void updateMappedSize(Mat_<apap_float> H, double x, double y,
        Point2d& corner, Point2d& br);

    /* Parameters */
    int cell_height_;
    int cell_width_;
    double gamma_;
    double sigma_;

    /* Temp variables */
    Mat_<apap_float> A_;
    Mat_<apap_float> W_;

    Mat_<apap_float>** H_s_;
    Mat_<apap_float>** H_r_;
    int cell_rows_, cell_cols_;

    Mat Hg_;
    Mat S_;

    Point2d K_min_, K_max_, K_1_, K_2_;

    vector<Point2d> anchor_points;
    Point2d offset;

    Mat T_src_;
    Mat T_dst_;


};

#endif //AANAP_APAP_H