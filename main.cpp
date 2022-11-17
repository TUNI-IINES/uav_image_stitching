#include <opencv2/opencv.hpp>
#include <vector>
#include <iostream>
#include "pipeline.cpp"

using namespace cv;
using namespace std;

int main(int argc, char** argv) {
    if (argc > 0) {

    	std::vector<cv::Mat> src;
        std::vector<double> time_prof1 {0.0, 0.0, 0.0, 0.0};
        std::vector<double> time_prof2 {0.0, 0.0, 0.0, 0.0};
        
        // With vscode, the target dirs is the build folder, so change this as ../imgs/
        string fov_0_path = "../generated_img_300/fov_0_" + string(argv[1]) + ".jpg"; 
        string fov_1_path = "../generated_img_300/fov_1_" + string(argv[1]) + ".jpg"; 
        string fov_2_path = "../generated_img_300/fov_2_" + string(argv[1]) + ".jpg"; 
        string fout_path = "../stitch_output/stitch_" + string(argv[1]) + ".jpg";

        cout << "Processing for " << string(argv[1]) << "ms" << endl;
    
        Mat fov_0 = imread(fov_0_path);
        Mat fov_1 = imread(fov_1_path);
        Mat fov_2 = imread(fov_2_path);
    
        // Stitch fov 0 and 1
    	src.push_back(fov_0);
    	src.push_back(fov_1);
        Mat stitched_image = pipeline(src, time_prof1);
        src.clear();
    
        // Stitch again with fov 2
    	src.push_back(stitched_image);
        src.push_back(fov_2);
        stitched_image = pipeline(src, time_prof2);
    
        printf("1st stitch in seconds: %.3f %.3f %.3f %.3f\n", time_prof1[0], time_prof1[1], time_prof1[2], time_prof1[3]);
        printf("2nd stitch in seconds: %.3f %.3f %.3f %.3f\n", time_prof2[0], time_prof2[1], time_prof2[2], time_prof2[3]);
        double total = 0.0;
        for (auto& n : time_prof1) total += n;
        for (auto& n : time_prof2) total += n;
        printf("total processing time in seconds: %.3f\n", total);

    	imwrite(fout_path, stitched_image); // Save image before further video resizing
        
        // std::ofstream myfile;
        // myfile.open ("time_profiling.csv");
        std::ofstream myfile( "time_profiling.csv", std::ios::app ) ;
        myfile << argv[1] << ", ,";
        myfile << time_prof1[0] << "," << time_prof1[1] << "," << time_prof1[2] << "," << time_prof1[3]<< ", ,";
        myfile << time_prof2[0] << "," << time_prof2[1] << "," << time_prof2[2] << "," << time_prof2[3]<< ", ,";
        myfile << total << "\n";
        myfile.close();
    }

    return 0;
}