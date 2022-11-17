#include <opencv2/opencv.hpp>
#include <vector>
#include <iostream>
#include "ImageSimulator.cpp"
#include "param.h"

using namespace cv;
using namespace std;

int main(int argc, char** argv) {
    // Read an image
    string image_path = "../imgs/b1.JPG"; // With vscode, the target dirs is the build folder, so change this as ../imgs/
    Mat canvas_img = imread(image_path);
    
    int height = canvas_img.size[0];
    int width = canvas_img.size[1];

    // Given a position and a radius, compute a rectangle with a specified image ratio. Then crop parts of the image
    // Allow multiple points & radius (for multiple camera views)
    // vector<Point2d> points = { Point2d(800, 800), Point2d(1300, 800), Point2d(1800, 800) };
    // vector<float> radius = { 500.0f, 500.0f, 500.0f };
    int idx_t = 0;
    vector<Point2d> points = { 
        Point2d( OR_PXL_X + floor(MTR2PXL*TRAJECTORY[idx_t][1]), OR_PXL_Y - floor(MTR2PXL*TRAJECTORY[idx_t][2])  ), 
        Point2d( OR_PXL_X + floor(MTR2PXL*TRAJECTORY[idx_t][5]), OR_PXL_Y - floor(MTR2PXL*TRAJECTORY[idx_t][6])  ), 
        Point2d( OR_PXL_X + floor(MTR2PXL*TRAJECTORY[idx_t][9]), OR_PXL_Y - floor(MTR2PXL*TRAJECTORY[idx_t][10]) )
    };
    vector<float> radius = { float(MTR2PXL*(TRAJECTORY[idx_t][4])), float(MTR2PXL*(TRAJECTORY[idx_t][8])), float(MTR2PXL*(TRAJECTORY[idx_t][12])) };

    // Add points and radius into Circle
    vector<Circle*> circles;
    for (int i = 0; i < points.size(); i++) {
        Circle* new_circle = new Circle(i, points[i], radius[i]);
        circles.push_back(new_circle);
    }

    int it = 0, cn = 0;
    int col_order[3] = { 0, 4, 8};
    int pos_x = 0, pos_y = 0;
    float temp_rad = -1.0f, rad = -1.0f;
    int time_ms = 0;

    ImageSimulator image_simulator = ImageSimulator(canvas_img);
    // image_simulator.load_frames();
    image_simulator.init_frames(SIMULATOR_FPS);

	// Iterate over time 
	for (int i = 9; i < SIMULATOR_MAX_TIME; i++) {
        // Iterate over fps counter
    	for (int fps = 0; fps < SIMULATOR_FPS; fps++) {
            it = i*50 + fps;
            time_ms = int(TRAJECTORY[it][0] * 1000);

            cout << i << ' ' << fps << ' ' << it << ' ' << time_ms << endl;

            cn = 0;
    	    for (Circle* circle : circles) {
                pos_x = OR_PXL_X + floor(MTR2PXL*TRAJECTORY[it][col_order[cn]+1]); 
                pos_y = OR_PXL_Y - floor(MTR2PXL*TRAJECTORY[it][col_order[cn]+2]);
                rad = (MTR2PXL*(TRAJECTORY[it][col_order[cn]+4]));
    	    	// Update the circles first
                circle->update_pose( pos_x, pos_y, rad);
                cn += 1;
            }
            image_simulator.stitch_frames(circles, time_ms);
        }
    }

    image_simulator.close_simulator();

    return 0;
}