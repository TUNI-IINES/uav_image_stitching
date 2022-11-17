#include "pipeline.h"

#define IMG_NUM 2

vector<Mat> xmaps_, ymaps_, final_warped_masks_;
Rect dst_roi_;
vector<Point> corners_;
vector<Size> sizes_;
const int parallel_num_ = 4;
vector<Mat> final_blend_masks_;

vector<Mat> blend_weight_maps_;

void meshgrid(const cv::Range& xgv, const cv::Range& ygv, cv::Mat& X, cv::Mat& Y)
{
    std::vector<int> t_x, t_y;
    for (int i = xgv.start; i <= xgv.end; i++) t_x.push_back(i);
    for (int j = ygv.start; j <= ygv.end; j++) t_y.push_back(j);

    cv::repeat(cv::Mat(t_x).t(), int(t_y.size()), 1, X);
    cv::repeat(cv::Mat(t_y), 1, int(t_x.size()), Y);
}

void circshift(Mat& out, const Point& delta)
{
    Size sz = out.size();

    // error checking
    assert(sz.height > 0 && sz.width > 0);

    // no need to shift
    if ((sz.height == 1 && sz.width == 1) || (delta.x == 0 && delta.y == 0))
        return;

    // delta transform
    int x = delta.x;
    int y = delta.y;
    if (x > 0) x = x % sz.width;
    if (y > 0) y = y % sz.height;
    if (x < 0) x = x % sz.width + sz.width;
    if (y < 0) y = y % sz.height + sz.height;


    // in case of multiple dimensions
    vector<Mat> planes;
    split(out, planes);

    for (size_t i = 0; i < planes.size(); i++)
    {
        // vertical
        Mat tmp0, tmp1, tmp2, tmp3;
        Mat q0(planes[i], Rect(0, 0, sz.width, sz.height - y));
        Mat q1(planes[i], Rect(0, sz.height - y, sz.width, y));
        q0.copyTo(tmp0);
        q1.copyTo(tmp1);
        tmp0.copyTo(planes[i](Rect(0, y, sz.width, sz.height - y)));
        tmp1.copyTo(planes[i](Rect(0, 0, sz.width, y)));

        // horizontal
        Mat q2(planes[i], Rect(0, 0, sz.width - x, sz.height));
        Mat q3(planes[i], Rect(sz.width - x, 0, x, sz.height));
        q2.copyTo(tmp2);
        q3.copyTo(tmp3);
        tmp2.copyTo(planes[i](Rect(x, 0, sz.width - x, sz.height)));
        tmp3.copyTo(planes[i](Rect(0, 0, x, sz.height)));
    }

    merge(planes, out);
}

void fftshift(Mat& out)
{
    Size sz = out.size();
    Point pt(0, 0);
    pt.x = (int)floor(sz.width / 2.0);
    pt.y = (int)floor(sz.height / 2.0);
    circshift(out, pt);
}

void ifft2(const Mat& src, Mat& Fourier)
{
    int mat_type = src.type();
    assert(mat_type < 15); //Unsupported Mat datatype

    if (mat_type < 7)
    {
        Mat planes[] = { Mat_<double>(src), Mat::zeros(src.size(),CV_64F) };
        merge(planes, 2, Fourier);
        dft(Fourier, Fourier, DFT_INVERSE + DFT_SCALE, 0);
    }
    else // 7<mat_type<15
    {
        Mat tmp;
        dft(src, tmp, DFT_INVERSE + DFT_SCALE, 0);
        //vector<Mat> planes;
        //split(tmp, planes);

        //magnitude(planes[0], planes[1], planes[0]); //Change complex to magnitude
        //Fourier = planes[0];
        Fourier = tmp;
    }
}

void fft2(const Mat& src, Mat& Fourier)
{
    int mat_type = src.type();
    assert(mat_type < 15); //Unsupported Mat datatype

    if (mat_type < 7)
    {
        Mat planes[] = { Mat_<double>(src), Mat::zeros(src.size(),CV_64F) };
        merge(planes, 2, Fourier);
        dft(Fourier, Fourier);
    }
    else // 7<mat_type<15
    {
        Mat tmp;
        dft(src, tmp);
        vector<Mat> planes;
        split(tmp, planes);
        magnitude(planes[0], planes[1], planes[0]); //Change complex to magnitude
        Fourier = planes[0];
    }
}

Mat createGabor(vector<int> orientationPerScale, int img_rows_w_padding, int img_cols_w_padding) {
    double pi = 3.14159265358979;
    int N_scales = int(orientationPerScale.size());
    int N_filters = orientationPerScale[0] * N_scales; // if all elements in orientationPerScale are the same

    int l = 0;
    vector<vector<double>> params;

    for (int i = 0; i < N_scales; i++) {
        for (int j = 0; j < orientationPerScale[i]; j++) {
            l += 1;
            params.push_back({0.35, 0.3/(pow(1.85, (i))), 16.0*pow(double(orientationPerScale[i]), 2)/pow(32, 2), (pi/orientationPerScale[i])*(j)});
        }
    }

    // Frequencies
    Mat fx, fy;
    meshgrid(Range(int(ceil(-img_cols_w_padding / 2.0)), int(ceil(img_cols_w_padding / 2.0 - 1))), Range(int(ceil(-img_rows_w_padding / 2.0)), int(ceil(img_rows_w_padding / 2.0 - 1))), fx, fy);
    fx.convertTo(fx, CV_64F);
    fy.convertTo(fy, CV_64F);

    // fftshift
    Mat fr, t, fx_pow, fy_pow;
    /*cv::multiply(fx, fx, fx_pow);
    cv::multiply(fy, fy, fy_pow);*/
    cv::pow(fx, 2, fx_pow);
    cv::pow(fy, 2, fy_pow);
    //fr.convertTo(fr, CV_64F);
    cv::sqrt(fx_pow + fy_pow, fr);
    fftshift(fr);

    //t = fx + sqrt(-1) * fy;
    //t.convertTo(t, CV_64F);
    cv::phase(fx, fy, t);
    cv::subtract(t, 2 * pi, t, (t > pi));
    fftshift(t);

    // Transfer function
    Mat G = Mat::zeros(img_rows_w_padding, img_cols_w_padding, CV_64FC(N_filters));
    vector<Mat> G_channels;
    split(G, G_channels);
    for (int i = 0; i < N_filters; i++) {
        Mat tr = t + params[i][3];
        Mat op1 = tr < -pi;
        op1.convertTo(op1, CV_64F, 1.0/255.0); // Because op1 and op2 are 0 and 255.0
        Mat op2 = tr > pi;
        op2.convertTo(op2, CV_64F, 1.0/255.0);
        tr = tr + 2 * pi * op1 - 2 * pi * op2;

        Mat B, C;
        Mat A = fr / img_cols_w_padding / params[i][1] - 1;
        cv::pow(A, 2, B);
        cv::pow(tr, 2, C);
        Mat result = -10 * params[i][0] * B - 2 * params[i][2] * pi * C;
        cv::exp(result, G_channels[i]);
    }
    merge(G_channels, G);

    return G;
    //cout << "Done" << endl;
}

Mat padArray(Mat img, int rowPad, int colPad) {
    // This is "symmetric" case in the Matlab's padArray (other option is different)
    int n = img.rows;
    int m = img.cols;
    Mat temp1 = Mat::zeros(n, m + colPad * 2, img.type());
    Mat temp2 = Mat::zeros(n + rowPad * 2, m + colPad * 2, img.type());
    for (int i = 0; i < colPad; i++) {
        img.col(i).copyTo(temp1.col(colPad - 1 - i));
        img.col(m - 1 - i).copyTo(temp1.col(m + colPad + i));
    }
    img.copyTo(temp1.colRange(colPad, m + colPad));

    for (int j = 0; j < rowPad; j++) {
        temp1.row(j).copyTo(temp2.row(rowPad - 1 - j));
        temp1.row(n - 1 - j).copyTo(temp2.row(n + rowPad + j));
    }
    temp1.copyTo(temp2.rowRange(rowPad, n + rowPad));
    return temp2;
}

Mat GaborConvole(Mat gray_img, Mat gabor_filter, int boundary_extension) {
    int ny = gabor_filter.rows;
    int nx = gabor_filter.cols;
    int N_filters = gabor_filter.channels();
    int img_row = gray_img.rows;
    int img_col = gray_img.cols;
    // Pad image
    Mat pad_gray_img = padArray(gray_img, boundary_extension, boundary_extension);
    Mat res;
    fft2(pad_gray_img, res);
    //dft(pad_gray_img, res, DFT_SCALE | DFT_COMPLEX_OUTPUT);
    // Convolving
    vector<Mat> gabor_filter_vect;
    split(gabor_filter, gabor_filter_vect);
    
    Mat texture_features = Mat::zeros(img_row, img_col, CV_64FC(N_filters));
    vector<Mat> texture_features_vect;
    split(texture_features, texture_features_vect);

    for (int i = 0; i < N_filters; i++) {
        Mat rep_mat = repeat(gabor_filter_vect[i], 1, 1);
        // Convert rep_mat to a 3-channel Mat
        Mat rep_mat_fc2;
        Mat t[] = {rep_mat, rep_mat };
        merge(t, 2, rep_mat_fc2);
        Mat ig;
        cv::multiply(res, rep_mat_fc2, ig);
        //multiply(res, rep_mat, ig);
        Mat inverse_transform;
        //dft(ig, inverse_transform, DFT_INVERSE);
        ifft2(ig, inverse_transform);
        vector<Mat> inverse_transform_vect;
        split(inverse_transform, inverse_transform_vect);
        Mat mag(inverse_transform.size(), CV_64F);
        magnitude(inverse_transform_vect[0], inverse_transform_vect[1], mag);
        
        Mat b = mag(Range(boundary_extension, ny - boundary_extension), Range(boundary_extension, nx - boundary_extension));
        texture_features_vect[i] = b;
    }
    merge(texture_features_vect, texture_features);
    return texture_features;
}

Mat getGaborResponse(Mat input_gray_img) {
    int norient = 8; // Number of wavelet scales
    int nscale = 5; // Number of filter orientations
    int boundary_extension = 32; // Number of pixels to pad

    vector<int> orientationPerScale(nscale, norient);
    int img_rows = input_gray_img.rows;
    int img_cols = input_gray_img.cols;

    // Matlab: GaborFilters = createGabor(orientationsPerScale, [ImgRow ImgCol] +2*boundaryExtension);
    Mat G = createGabor(orientationPerScale, img_rows + 2 * boundary_extension, img_cols + 2 * boundary_extension);

    Mat texture_features = GaborConvole(input_gray_img, G, boundary_extension);
    return texture_features;
}

class IndexCompare {
private:
    Mat _data;
public:
    IndexCompare(Mat data) :_data(data)
    {}

    bool operator()(int i, int j) const
    {
        return _data.at<double>(i, 0) < _data.at<double>(j, 0);
    }
};


pair<Mat, Mat> top(Mat X, int k) {
    int d = X.rows;
    int n = X.cols;

    Mat values = Mat(k, n, CV_64F);
    Mat indices = Mat(k, n, CV_8U);
    
    vector<int> tmp(d);
    
    for (int i = 0; i < n; i++) {
        
        for (int j = 0; j < d; j++)
            tmp[j] = j;

        partial_sort(tmp.begin(), tmp.begin() + k, tmp.end(), IndexCompare(X.col(i)));
        
        for (int j = 0; j < k; j++) {
            values.at<double>(j, i) = X.at<double>(tmp[j], i);
            indices.at<uchar>(j, i) = tmp[j];
        }
    }
    pair<Mat, Mat> result = make_pair(values, indices);
    return result;
}

pair<Mat, Mat> knn(Mat img_texture_response, Mat texton, int k) {
    // Perform element-wise plus arrays
    Mat A = (-2) * texton.t() * img_texture_response;
    vector<double> result_vect;
    for (int i = 0; i < texton.cols; i++) {
        double result = texton.col(i).dot(texton.col(i));
        result_vect.push_back(result);
    }
    Mat B = Mat(result_vect);
    for (int i = 0; i < A.cols; i++) {
        // Plus col-by-col
        A.col(i) = A.col(i) + B;
    }
    pair<Mat, Mat> knn_result = top(A, k);
    Mat distance = knn_result.first;
    Mat index = knn_result.second;

    vector<double> result_vect_1;
    for (int i = 0; i < img_texture_response.cols; i++) {
        double result = img_texture_response.col(i).dot(img_texture_response.col(i));
        result_vect_1.push_back(result);
    }
    Mat B_1 = Mat(result_vect_1);
    transpose(B_1, B_1);
    distance = distance + B_1;

    pair<Mat, Mat> result = make_pair(distance, index);

    return result;
}

Mat getInputTextureImg(Mat texture_features) {
    ifstream texton_file("texton.txt");
    int texton_rows = 64;
    int texton_cols = 40;
    Mat texton = Mat::zeros(texton_rows, texton_cols, CV_64F);
    for (int i = 0; i < texton_rows; i++) {
        for (int j = 0; j < texton_cols; j++) {
            texton_file >> texton.at<double>(i, j);
        }
    }
    int img_row = texture_features.rows;
    int img_col = texture_features.cols;
    int num_dim = texture_features.channels();

    int num_texton = texton.rows;
    int num_img_pixels = img_row * img_col;
    // Reshape texture_features
    Mat img_texture_response = texture_features.reshape(1, num_img_pixels);
    // Transpose
    transpose(img_texture_response, img_texture_response);
    transpose(texton, texton);
    // Perform knn
    pair<Mat, Mat> knn_result = knn(img_texture_response, texton, 1);
    Mat fast_dist_matrix = knn_result.first;
    Mat fast_min_idx = knn_result.second;

    Mat tmp_fast_min_idx = fast_min_idx.row(0);
    Mat fast_texton_img = tmp_fast_min_idx.reshape(1, img_row);

    return fast_texton_img;
}

Mat getTextureFeatures(Mat input_gray_img) {
    Mat gray_img;
    // Smoothing before getting response
    GaussianBlur(input_gray_img, gray_img, Size(5, 5), 1.0);

    // Doing for the input image
    Mat input_texture_features = getGaborResponse(gray_img);
    // Getting the texture image based on texton
    //Mat input_texture_img = getInputTextureImg(input_texture_features);

    return input_texture_features;
}

Mat Dx(Mat u) {
    int rows = u.rows;
    int cols = u.cols;
    int p = u.channels();
    Mat d;

    if (p == 1) {
        d.create(u.size(), CV_64F);
    }
    else if (p == 3) {
        d.create(u.size(), CV_64FC3);
    }
    d.setTo(Scalar::all(0.0));

    /*for (int j = 1; j < cols; j++) {
        for (int i = 0; i < rows; i++) {
            d.at<double>(i, j) = u.at<double>(i, j) - u.at<double>(i, j - 1);
        }
    }*/
    d(Range(0, rows), Range(1, cols)) = u(Range(0, rows), Range(1, cols)) - u(Range(0, rows), Range(0, cols - 1));
    d.col(0) = u.col(0) - u.col(cols-1);
    /*for (int i = 0; i < rows; i++) {
        d.at<double>(i, 0) = u.at<double>(i, 0) - u.at<double>(i, cols - 1);
    }*/
    return d;
}

Mat Dy(Mat u) {
    int rows = u.rows;
    int cols = u.cols;
    int p = u.channels();
    Mat d;

    if (p == 1) {
        d.create(u.size(), CV_64F);
    }
    else if (p == 3) {
        d.create(u.size(), CV_64FC3);
    }
    d.setTo(Scalar::all(0.0));

    /*for (int i = 1; i < rows; i++) {

        for (int j = 0; j < cols; j++) {
            d.at<double>(i, j) = u.at<double>(i, j) - u.at<double>(i - 1, j);
        }
    }
    for (int j = 0; j < cols; j++) {
        d.at<double>(0, j) = u.at<double>(0, j) - u.at<double>(rows - 1, j);
    }*/
    d(Range(1, rows), Range(0, cols)) = u(Range(1, rows), Range(0, cols)) - u(Range(0, rows - 1), Range(0, cols));
    d.row(0) = u.row(0) - u.row(rows - 1);
    return d;
}


pair<Mat, Mat> cost(Mat labels, Mat src_img, Mat dst_img, int sp_num) {
    vector<Mat> src_channels, dst_channels;
    split(src_img, src_channels);
    split(dst_img, dst_channels);
    Mat Y1 = src_channels[0];
    Mat U1 = src_channels[1];
    Mat V1 = src_channels[2];
    Mat Y2 = dst_channels[0];
    Mat U2 = dst_channels[1];
    Mat V2 = dst_channels[2];

    int m = src_img.size[0];
    int n = src_img.size[1];

    double w_y = 0.9;
    double w_u = 0.05;
    double w_v = 0.05;
    double lambda_c = 1.0;
    double lambda_g = 1.5;

    Mat cost_intensity(src_img.size(), CV_64F, Scalar(0));
    Mat total_cost(Size(sp_num, sp_num), CV_64F, Scalar(0));
    //Mat cost_in(Size(sp_num, sp_num), CV_64F, Scalar(0));
    //Mat cost_gr(Size(sp_num, sp_num), CV_64F, Scalar(0));
    //Mat cost_ig(Size(sp_num, sp_num), CV_64F, Scalar(0));
    //Mat sp_intensity(Size(sp_num, 1), CV_64F, Scalar(0));
    //Mat sp_gradient(Size(sp_num, 1), CV_64F, Scalar(0));
    //Mat sp_ig(Size(sp_num, 1), CV_64F, Scalar(0));
    //Mat sp_cost(Size(sp_num, 1), CV_64F, Scalar(0));

    // Compute the adjacency matrix
    Mat adjacency(Size(sp_num, sp_num), CV_8U, Scalar(0));
    for (int i = 1; i < m - 1; i++) {
        int* labels_ptr_current = labels.ptr<int>(i);
        int* labels_ptr_next = labels.ptr<int>(i+1);
        int* labels_ptr_prev = labels.ptr<int>(i-1);
        for (int j = 1; j < n - 1; j++) {
            //cout << "row: " + to_string(i) + ", col: " + to_string(j) << endl;
            /*if (i == 328 && j == 448) {
                cout << "debug here" << endl;
            }*/
            if (labels_ptr_current[j] != labels_ptr_current[j+1] && labels_ptr_current[j] != 0 && labels_ptr_current[j + 1] != 0) {
                adjacency.at<uchar>(labels_ptr_current[j] - 1, labels_ptr_current[j + 1] - 1) = 1;
                adjacency.at<uchar>(labels_ptr_current[j + 1] - 1, labels_ptr_current[j] - 1) = 1;
            }
            if (labels_ptr_current[j] != labels_ptr_current[j-1] && labels_ptr_current[j] != 0 && labels_ptr_current[j - 1] != 0) {
                adjacency.at<uchar>(labels_ptr_current[j] - 1, labels_ptr_current[j - 1] - 1) = 1;
                adjacency.at<uchar>(labels_ptr_current[j - 1] - 1, labels_ptr_current[j] - 1) = 1;
            }
            if (labels_ptr_current[j] != labels_ptr_next[j] && labels_ptr_current[j] != 0 && labels_ptr_next[j] != 0) {
                adjacency.at<uchar>(labels_ptr_current[j] - 1, labels_ptr_next[j] - 1) = 1;
                adjacency.at<uchar>(labels_ptr_next[j] - 1, labels_ptr_current[j] - 1) = 1;
            }
            if (labels_ptr_current[j] != labels_ptr_prev[j] && labels_ptr_current[j] != 0 && labels_ptr_prev[j] != 0) {
                adjacency.at<uchar>(labels_ptr_current[j] - 1, labels_ptr_prev[j] - 1) = 1;
                adjacency.at<uchar>(labels_ptr_prev[j] - 1, labels_ptr_current[j] - 1) = 1;
            }
            if (labels_ptr_current[j] != labels_ptr_next[j+1] && labels_ptr_current[j] != 0 && labels_ptr_next[j + 1] != 0) {
                adjacency.at<uchar>(labels_ptr_current[j] - 1, labels_ptr_next[j + 1] - 1) = 1;
                adjacency.at<uchar>(labels_ptr_next[j + 1] - 1, labels_ptr_current[j] - 1) = 1;
            }
            if (labels_ptr_current[j] != labels_ptr_next[j-1] && labels_ptr_current[j] != 0 && labels_ptr_next[j - 1] != 0) {
                adjacency.at<uchar>(labels_ptr_current[j] - 1, labels_ptr_next[j - 1] - 1) = 1;
                adjacency.at<uchar>(labels_ptr_next[j - 1] - 1, labels_ptr_current[j] - 1) = 1;
            }
            if (labels_ptr_current[j] != labels_ptr_prev[j+1] && labels_ptr_current[j] != 0 && labels_ptr_prev[j + 1] != 0) {
                adjacency.at<uchar>(labels_ptr_current[j] - 1, labels_ptr_prev[j + 1] - 1) = 1;
                adjacency.at<uchar>(labels_ptr_prev[j + 1] - 1, labels_ptr_current[j] - 1) = 1;
            }
            if (labels_ptr_current[j] != labels_ptr_prev[j-1] && labels_ptr_current[j] != 0 && labels_ptr_prev[j - 1] != 0) {
                adjacency.at<uchar>(labels_ptr_current[j] - 1, labels_ptr_prev[j - 1] - 1) = 1;
                adjacency.at<uchar>(labels_ptr_prev[j - 1] - 1, labels_ptr_current[j] - 1) = 1;
            }
        }
    }

    //for (int i = 0; i < m; i++) {
    //    double* cost_intensity_ptr = cost_intensity.ptr<double>(i);
    //    double* Y1_ptr = Y1.ptr<double>(i);
    //    double* Y2_ptr = Y2.ptr<double>(i);
    //    double* U2_ptr = U2.ptr<double>(i);
    //    double* V2_ptr = V2.ptr<double>(i);
    //    double* U1_ptr = U1.ptr<double>(i);
    //    double* V1_ptr = V1.ptr<double>(i);

    //    for (int j = 0; j < n; j++) {
    //        //cout << "row: " + to_string(i) + ", col: " + to_string(j) << endl;
    //        cost_intensity_ptr[j] = w_y * abs(Y1_ptr[j] - Y2_ptr[j]) 
    //            + w_u * abs(U1_ptr[j] - U2_ptr[j]) 
    //            + w_v * abs(V1_ptr[j] - V2_ptr[j]);
    //    }
    //}

    cost_intensity = w_y * cv::abs(Y1 - Y2) + w_u * cv::abs(U1 - U2) + w_v * cv::abs(V1 - V2);
    
    Mat cost_gradient, cost_texture;
    cost_gradient = abs(Dx(Y1) - Dx(Y2)) + abs(Dy(Y1) - Dy(Y2));
    cost_texture = getTextureFeatures(Y1) + getTextureFeatures(Y2);
    cost_intensity = cost_intensity.reshape(cost_intensity.channels(), m * n);
    cost_gradient = cost_gradient.reshape(cost_gradient.channels(), m * n);
    cost_texture = cost_texture.reshape(1, m * n);

    //double minVal_1;
    //double maxVal_1;
    //Point minLoc_1;
    //Point maxLoc_1;

    //minMaxLoc(cost_texture, &minVal_1, &maxVal_1, &minLoc_1, &maxLoc_1);

    //cout << "min val: " << minVal_1 << endl;
    //cout << "max val: " << maxVal_1 << endl;

    Mat cost_texture_n = Mat(cost_texture.rows, 1, cost_texture.type());
    for (int i = 0; i < cost_texture.rows; i++) {
        cost_texture_n.at<double>(i, 0) = norm(cost_texture.row(i), NORM_L2);
    }
    // Normalize of texture
    Mat cost_texture_nor;
    cv::normalize(cost_texture_n, cost_texture_nor, 1.0, 0.0, NORM_MINMAX);
    Mat nor_texture = cost_texture_nor.reshape(cost_texture_nor.channels(), m);

    Mat p_cost;
    cv::multiply((lambda_c * cost_intensity + lambda_g * cost_gradient), cost_texture_nor, p_cost);
    Mat p_r_cost = p_cost.reshape(p_cost.channels(), m);

    vector<double> sp_cost;
    for (int i = 0; i < sp_num; i++) {
        Mat sp = (labels == i+1);
        Mat currSP_idx = (sp == 255);

        // Compute the mean feature of each superpixel
        Scalar cost = mean(p_r_cost, currSP_idx).val[0];
        sp_cost.push_back(cost.val[0]);
    }

    for (int i = 0; i < sp_num; i++) {
        uchar* adjacency_ptr = adjacency.ptr<uchar>(i);
        double* total_cost_ptr = total_cost.ptr<double>(i);
        for (int j = 0; j < sp_num; j++) {
            if (adjacency_ptr[j] == 1) {
                total_cost_ptr[j] = sp_cost[i] + sp_cost[j];
            }
        }
    }
    pair<Mat, Mat> final_result;
    final_result = make_pair(total_cost, p_r_cost);
    return final_result;
}

Mat isEdge(Mat image) {
    Size s = image.size();
    int rows = image.rows;
    int cols = image.cols;

    Mat edge = Mat::zeros(image.size(), CV_8U);

    for (int i = 1; i < rows-1; i++) {
        uchar *image_ptr_current = image.ptr<uchar>(i);
        uchar* image_ptr_next = image.ptr<uchar>(i + 1);
        uchar* image_ptr_prev = image.ptr<uchar>(i - 1);
        uchar* edge_ptr = edge.ptr<uchar>(i);
        for (int j = 1; j < cols-1; j++) {
            //cout << "i: " + to_string(i) + " j: " + to_string(j) << endl;
            if (image_ptr_current[j] != image_ptr_next[j]) {
                edge_ptr[j] = 255;
            }
            else if (image_ptr_current[j] != image_ptr_prev[j]) {
                edge_ptr[j] = 255;
            }
            else if (image_ptr_current[j] != image_ptr_current[j+1]) {
                edge_ptr[j] = 255;
            }
            else if (image_ptr_current[j] != image_ptr_current[j-1]) {
                edge_ptr[j] = 255;
            }
        }
    }
    //imshow("Edge", edge);
    //waitKey(0);
    return edge;
}

template <typename T>
vector<size_t> sort_indexes(const vector<T>& v) {

    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values 
    stable_sort(idx.begin(), idx.end(),
        [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

    return idx;
}

Mat image_blending_color(vector<Mat> images_warped, vector<Mat> seamed_masks, int sp_num, Mat labels, Mat cutline, Mat labelclass) {
    int beta = 50;

    int m = images_warped[0].rows;
    int n = images_warped[0].cols;
    Mat img1_double, img2_double;
    images_warped[0].convertTo(img1_double, CV_64FC3);
    images_warped[1].convertTo(img2_double, CV_64FC3);
    vector<Mat> img1_channels, img2_channels;
    split(img1_double, img1_channels);
    split(img2_double, img2_channels);

    // Get the color difference of the seam
    int num = 1;
    vector<Vec3d> edge_1, edge_2, bias;
    for (int i = 0; i < m; i++) {
        Vec3d* img1_double_ptr = img1_double.ptr<Vec3d>(i);
        Vec3d* img2_double_ptr = img2_double.ptr<Vec3d>(i);
        for (int j = 0; j < n; j++) {
            if (cutline.at<int>(i, j) == 1) {
                edge_1.push_back(img1_double_ptr[j]);
                edge_2.push_back(img2_double_ptr[j]);
                bias.push_back(img1_double_ptr[j] - img2_double_ptr[j]);
                num += 1;
            }
        }
    }
    num -= 1;

    // Randomly select a point in the superpixel to represent the superpixel
    vector<Vec3d> candidate;
    Mat count = Mat::zeros(sp_num, 3, CV_64F);

    for (int k = 0; k < sp_num; k++) {
        count.at<double>(k, 0) = 1;
        map<double, int> tmp1, tmp2, tmp3;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                if (labels.at<int>(i, j) == k + 1 && seamed_masks[1].at<uchar>(i, j) == 255) {
                    /*tmp1.push_back(img2_channels[0].at<double>(i, j));
                    tmp2.push_back(img2_channels[1].at<double>(i, j));
                    tmp3.push_back(img2_channels[2].at<double>(i, j));*/
                    // tmp1
                    if (tmp1.find(img2_channels[0].at<double>(i, j)) == tmp1.end()) {
                        // not found
                        tmp1.insert(make_pair(img2_channels[0].at<double>(i, j), 1));
                    }
                    else {
                        tmp1.at(img2_channels[0].at<double>(i, j)) += 1;
                    }
                    // tmp2
                    if (tmp2.find(img2_channels[1].at<double>(i, j)) == tmp2.end()) {
                        // not found
                        tmp2.insert(make_pair(img2_channels[1].at<double>(i, j), 1));
                    }
                    else {
                        tmp2.at(img2_channels[1].at<double>(i, j)) += 1;
                    }
                    // tmp3
                    if (tmp3.find(img2_channels[2].at<double>(i, j)) == tmp3.end()) {
                        // not found
                        tmp3.insert(make_pair(img2_channels[2].at<double>(i, j), 1));
                    }
                    else {
                        tmp3.at(img2_channels[2].at<double>(i, j)) += 1;
                    }
                    count.at<double>(k, 0) += 1;
                }
            }
        }
        count.at<double>(k, 0) -= 1;
        std::map<double, int>::iterator max_tmp1 = max_element(
            tmp1.begin(), tmp1.end(), [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                return a.second < b.second;
            }
        );
        std::map<double, int>::iterator max_tmp2 = max_element(
            tmp2.begin(), tmp2.end(), [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                return a.second < b.second;
            }
        );
        std::map<double, int>::iterator max_tmp3 = max_element(
            tmp3.begin(), tmp3.end(), [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                return a.second < b.second;
            }
        );
        Vec3d maxes(max_tmp1->first, max_tmp2->first, max_tmp3->first);
        candidate.push_back(maxes);

    }
    // Calculate the relative weight
    //Mat color_R = Mat::zeros(sp_num, 3, CV_64F);
    vector<Vec3d> color_R(sp_num);
    
    for (int k = 0; k < sp_num; k++) {
        vector<double> w_r;
        for (int p = 0; p < num; p++) {
            Vec3d substraction = candidate[k] - edge_2[p];
            double delta = norm(substraction, NORM_L2);
            double w = exp(-pow(delta, 2) / (2 * pow(beta, 2)));
            w_r.push_back(w);
        }
        // Sort descendant
        vector<double> w_n1;
        vector<int> I_1;
        for (int idx : sort_indexes(w_r)) {
            w_n1.push_back(w_r[idx]);
            I_1.push_back(idx);
        }
        double num1 = round(num * 0.95);
        vector<double> w_n2(w_n1.begin(), w_n1.begin() + num1);
        double sum_w_n2 = 0.0;
        for (double w : w_n2) {
            sum_w_n2 += w;
        }
        if (sum_w_n2 != 0) {
            vector<double> w_n;
            for (double w : w_n2) {
                w_n.push_back(w / sum_w_n2);
            }
            for (int z = 0; z < num1; z++) {
                /*color_R.at<double>(k, 0) += w_n[z] * bias[I_1[z]](0);
                color_R.at<double>(k, 1) += w_n[z] * bias[I_1[z]](1);
                color_R.at<double>(k, 2) += w_n[z] * bias[I_1[z]](2);*/
                color_R[k] = color_R[k] + w_n[z] * bias[I_1[z]];
            }
        }
    }

    // Calculate the final color
    Mat pano = Mat::zeros(m, n, CV_64FC3);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            //cout << to_string(i) + " and " + to_string(j) << endl;
            /*if(i == 302 && j == 462) {
                cout << "debug" << endl;
            }*/
            if (labelclass.at<int>(i, j) == 1 || (seamed_masks[0].at<uchar>(i, j) == 255 && seamed_masks[1].at<uchar>(i, j) == 0)) {
                pano.at<Vec3d>(i, j) = img1_double.at<Vec3d>(i, j);
            }
            else if (labelclass.at<int>(i, j) == 2 || (seamed_masks[1].at<uchar>(i, j) == 255 && seamed_masks[0].at<uchar>(i, j) == 0)) {
                pano.at<Vec3d>(i, j) = img2_double.at<Vec3d>(i, j) + color_R[labels.at<int>(i, j)-1];
            }
            else {
                pano.at<Vec3d>(i, j) = { 255, 255, 255};
            }
        }
    }
    Mat pano_8U;
    pano.convertTo(pano_8U, CV_8UC3);
    //imshow("Pano", pano_8U);
    //waitKey(0);

    return pano_8U;
}

void spdisplay(Mat img, Mat cutline) {
    int m = img.rows;
    int n = img.cols;
    Mat output_img = img.clone();
    for (int i = 1; i < m - 1; i++) {
        Vec3b* img_ptr = output_img.ptr<Vec3b>(i);
        int* cutline_ptr = cutline.ptr<int>(i);
        for (int j = 1; j < n - 1; j++) {
            if (cutline_ptr[j] == 1) {
                img_ptr[j] = { 255, 0, 0 };
            }
        }
    }
    imshow("Image with cutline", output_img);
    waitKey(0);
}

void image_aligment(vector<Mat> &src, vector<Mat> &seamed_masks, vector<Mat> &images_warped, vector<Mat> &init_masks) {
    // Keypoint detection using SIFT
    Ptr<SIFT> sift = SIFT::create(0, 3, 0.01, 40);
    vector<ImageFeatures> features(IMG_NUM);
    for (int i = 0; i < IMG_NUM; ++i) {
        features[i].img_idx = i;
        features[i].img_size = src[i].size();
        sift->detect(src[i], features[i].keypoints);
        sift->compute(src[i], features[i].keypoints, features[i].descriptors);
    }

    vector<MatchesInfo> matches_info;
    BestOf2NearestMatcher matcher(false, 0.5f);
    matcher(features, matches_info);
    matcher.collectGarbage();

    corners_.resize(IMG_NUM);
    sizes_.resize(IMG_NUM);
    xmaps_.resize(IMG_NUM);
    ymaps_.resize(IMG_NUM);
    Ptr<AANAPWarper> warper = new AANAPWarper();
    warper->buildMaps(src, features, matches_info, xmaps_, ymaps_, corners_);

    // Get the mask of the overlapping area
    for (int i = 0; i < IMG_NUM; ++i)
        sizes_[i] = xmaps_[i].size();
    dst_roi_ = resultRoi(corners_, sizes_);
    final_warped_masks_.resize(IMG_NUM);

    for (int i = 0; i < IMG_NUM; i++)
    {
        init_masks[i].create(src[i].size(), CV_8U);
        init_masks[i].setTo(Scalar::all(255));
        remap(src[i], images_warped[i], xmaps_[i], ymaps_[i], INTER_LINEAR);
        remap(init_masks[i], final_warped_masks_[i], xmaps_[i], ymaps_[i], INTER_NEAREST, BORDER_CONSTANT);
        seamed_masks[i] = final_warped_masks_[i].clone();
        //imshow("seamed masks", images_warped[i]);
        //waitKey(0);
    }
}

void image_registration(Mat &mask, int m, int n, vector<Mat> &seamed_masks, vector<Mat> &images_warped, vector<Mat> &src_intensity, vector<Mat> &src_YUV, vector<Mat> &images_warped_clone, map<int, pair<int, int>> &intersection_list) {

    bitwise_and(seamed_masks[0], seamed_masks[1], mask);

    //imshow("Mask after getting overlapped area", mask);
    //waitKey(0);

    // Compute the intersection of images
    vector<Mat> edges;
    for (int i = 0; i < IMG_NUM; i++) {
        edges.push_back(isEdge(seamed_masks[i]));
    }

    Mat flag;
    bitwise_and(edges[0], edges[1], flag);
    //imshow("flag", flag);
    //waitKey(0);

    // Create intersection_list for the graph cut
    //map<int, pair<int, int>> intersection_list;
    for (int i = 0; i < m; i++) {
        uchar* flag_ptr = flag.ptr<uchar>(i);
        uchar* mask_ptr = mask.ptr<uchar>(i);
        for (int j = 0; j < n; j++) {
            if (flag_ptr[j] == 255 && mask_ptr[j] == 255) {
                intersection_list.insert({ i + j, make_pair(i, j) });
            }
        }
    }

    // Convert to YUV
    for (int i = 0; i < IMG_NUM; i++) {
        vector<Mat> channels(3);
        Mat Y, U, V;

        Mat img_intensity = Mat(images_warped[i].size(), CV_64F); // images_warped but has been scaled to 0-1
        Mat img_YUV = Mat(images_warped[i].size(), CV_64F);
        cv::normalize(images_warped[i], img_intensity, 0.0, 1.0, NORM_MINMAX, CV_64F);
        // Split the channels
        split(img_intensity, channels);
        // Access the channels
        Y = channels[2] * 0.299 + channels[1] * 0.587 + channels[0] * 0.114;
        U = channels[2] * -0.147 + channels[1] * -0.289 + channels[0] * 0.436;
        V = channels[2] * 0.615 + channels[1] * -0.515 + channels[0] * -0.100;
        // merge YUV channels
        vector<Mat> YUV_channels = { Y, U, V };
        merge(YUV_channels, img_YUV);
        // Push back the vectors
        src_intensity.push_back(img_intensity);
        src_YUV.push_back(img_YUV);
    }

    // Clone images_warped for SLIC algorithm
    for (Mat image : images_warped) {
        images_warped_clone.push_back(image.clone());
    }

    // Compute the valid area of each image
    for (int i = 0; i < m; i++) {
        uchar* mask_ptr = mask.ptr<uchar>(i);
        Vec3d* src_intensity_ptr = src_intensity[0].ptr<Vec3d>(i);
        Vec3b* image_warped_clone_ptr = images_warped_clone[0].ptr<Vec3b>(i);
        for (int j = 0; j < n; j++) {
            if (mask_ptr[j] != 255) {
                // Turn src_intensity[0]
                src_intensity_ptr[j] = { 0, 0, 0 };
                //// Turn images_warped[0]
                image_warped_clone_ptr[j] = { 0, 0, 0 };
            }
        }
    }
}

int slic_and_labels(int m, int n, int sp_size, float compactnessfactor, Mat& labels, Mat& images_warped_clone, Mat& mask, vector<Mat>& src_YUV, Mat& seamed_mask, bool pre) {
    // Smoothness cost

    // SLIC algorithm
    int sp_num = 0;

    Ptr<ximgproc::SuperpixelSLIC> slic = ximgproc::createSuperpixelSLIC(images_warped_clone, ximgproc::SLICO, sp_size, compactnessfactor);
    slic->iterate();
    slic->getLabels(labels);
    sp_num = slic->getNumberOfSuperpixels();

    labels = labels + 1;

    for (int i = 0; i < m; i++) {
        uchar* mask_ptr = mask.ptr<uchar>(i);
        uchar* seamed_mask_ptr = seamed_mask.ptr<uchar>(i);
        int* labels_ptr = labels.ptr<int>(i);
        for (int j = 0; j < n; j++) {
            if (pre) {
                if (mask_ptr[j] != 255) {
                    labels_ptr[j] = 0;
                    for (int k = 0; k < IMG_NUM; k++) {
                        src_YUV[k].at<Vec3d>(i, j) = { 0, 0, 0 };
                    }
                }
            }
            else {
                if (seamed_mask_ptr[j] != 255) {
                    labels_ptr[j] = 0;
                }
            }
        }
    }

    // Get valid area
    set<int> label_list;
    for (int i = 0; i < m; i++) {
        int* labels_ptr = labels.ptr<int>(i);
        for (int j = 0; j < n; j++) {
            label_list.insert(labels_ptr[j]);
        }
    }

    int label_idx = 0;
    for (auto element : label_list) {
        if (label_idx > 0) {
            for (int i = 0; i < m; i++) {
                int* labels_ptr = labels.ptr<int>(i);
                for (int j = 0; j < n; j++) {
                    if (labels_ptr[j] == element) {
                        labels_ptr[j] = label_idx;
                    }
                }
            }
        }
        label_idx++;
    }

    sp_num = int(label_list.size()) - 1;

    return sp_num;
}

Mat pipeline(vector<Mat> &src) {
    // Image aligment and registration

    vector<Mat> seamed_masks(IMG_NUM);
    vector<Mat> images_warped(IMG_NUM);
    vector<Mat> init_masks(IMG_NUM);

    image_aligment(src, seamed_masks, images_warped, init_masks);

    Mat mask = Mat::zeros(dst_roi_.size(), CV_8U);
    int m = dst_roi_.height;
    int n = dst_roi_.width;
    map<int, pair<int, int>> intersection_list;
    vector<Mat> src_intensity;
    vector<Mat> src_YUV;
    vector<Mat>images_warped_clone;

    image_registration(mask, m, n, seamed_masks, images_warped, src_intensity, src_YUV, images_warped_clone, intersection_list);

    // Smoothness cost

    Mat labels;
    int sp_num = slic_and_labels(m, n, 20, 5.0f, labels, images_warped_clone[0], mask, src_YUV, seamed_masks[0]);

    // Compute the cost matrix
    pair<Mat, Mat> costs = cost(labels, src_YUV[0], src_YUV[1], sp_num);
    Mat totalcost = costs.first;
    Mat p_r_cost = costs.second;
    // Matlab: Save `labels` and `totalcost` and move to C++.
    labels = labels - 1;

    // Graph cut

    int* mask_pointer = new int[mask.rows * mask.cols];
    int* label_ptr = new int[mask.rows * mask.cols];

    for (int i = 0; i < mask.rows; i++) {
        for (int j = 0; j < mask.cols; j++) {
            mask_pointer[i * mask.cols + j] = mask.at<uchar>(i, j)/255;
            label_ptr[i * mask.cols + j] = labels.at<int>(i, j); // Similar to when saving label.txt at Matlab
        }
    }

    Mat totalcost_float;
    totalcost.convertTo(totalcost_float, CV_32F);

    vector<float> totalcost_vect;
    if (totalcost_float.isContinuous()) {
        totalcost_vect.assign((float*)totalcost_float.data, (float*)totalcost_float.data + totalcost_float.total() * totalcost_float.channels());
    }
    else {
        for (int i = 0; i < totalcost_float.rows; ++i) {
            totalcost_vect.insert(totalcost_vect.end(), totalcost_float.ptr<float>(i), totalcost_float.ptr<float>(i) + totalcost_float.cols * totalcost_float.channels());
        }
    }
    // As can be seen from Matlab, the x, y coordinates at the intersection are reversed order
    pair<int, int> pts_1 = intersection_list.begin()->second;
    pair<int, int> pts_2 = intersection_list.rbegin()->second;
    
    pair<int*, int*> seamcut_result = seamcut(mask.cols, mask.rows, mask_pointer, label_ptr, sp_num, totalcost_vect, pts_1.second, pts_1.first, pts_2.second, pts_2.first);
    int* cutline_ptr = seamcut_result.first;
    int* labelclass_ptr = seamcut_result.second;

    // Transfer cutline_ptr and labelclass_ptr into cv::Mat
    Mat cutline(m, n, CV_32S, cutline_ptr);
    Mat labelclass(m, n, CV_32S, labelclass_ptr);

    // Post-process
    Mat labels_2;
    int sp_num_2 = slic_and_labels(m, n, 40, 7.0f, labels_2, images_warped_clone[1], mask, src_YUV, seamed_masks[1], false);

    // Image blending
    Mat pano = image_blending_color(images_warped, seamed_masks, sp_num_2, labels_2, cutline, labelclass);
    //spdisplay(pano, cutline);
    return pano;
}
