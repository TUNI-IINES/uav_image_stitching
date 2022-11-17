#include <iostream>
#include <algorithm>
#include <numeric>
#include <opencv2/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/xfeatures2d.hpp>
#include <opencv2/ximgproc/slic.hpp>
#include <math.h>
#include <chrono>

#include "aanap.cpp"
#include "seamcut.cpp"

void meshgrid(const Range& xgv, const Range& ygv, Mat& X, Mat& Y);
void circshift(Mat& out, const Point& delta);
void fftshift(Mat& out);
void ifft2(const Mat& src, Mat& Fourier);
void fft2(const Mat& src, Mat& Fourier);
Mat createGabor(vector<int> orientationPerScale, int img_rows_w_padding, int img_cols_w_padding);
Mat padArray(Mat img, int rowPad, int colPad);
Mat GaborConvole(Mat gray_img, Mat gabor_filter, int boundary_extension);
Mat getGaborResponse(Mat input_gray_img);
pair<Mat, Mat> top(Mat X, int k);
pair<Mat, Mat> knn(Mat img_texture_response, Mat texton, int k);
Mat getInputTextureImg(Mat texture_features);
Mat getTextureFeatures(Mat input_gray_img);
Mat Dx(Mat u);
Mat Dy(Mat u);
pair<Mat, Mat> cost(Mat labels, Mat src_img, Mat dst_img, int sp_num);
Mat isEdge(Mat image);
Mat image_blending_color(vector<Mat> images_warped, vector<Mat> seamed_masks, int sp_num, Mat labels, Mat cutline, Mat labelclass);
void spdisplay(Mat img, Mat cutline);
void image_aligment(vector<Mat>& src, vector<Mat>& seamed_masks, vector<Mat>& images_warped, vector<Mat>& init_masks);
void image_registration(Mat& mask, int m, int n, vector<Mat>& seamed_masks, vector<Mat>& images_warped, vector<Mat>& src_intensity, vector<Mat>& src_YUV, vector<Mat>& images_warped_clone, map<int, pair<int, int>>& intersection_list);
int slic_and_labels(int m, int n, int sp_size, float compactnessfactor, Mat& labels, Mat& images_warped_clone, Mat& mask, vector<Mat>& src_YUV, Mat& seamed_mask, bool pre = true);
Mat pipeline(vector<Mat> &src, std::vector<double> &time_prof);