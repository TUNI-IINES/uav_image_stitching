#include <opencv2/opencv.hpp>
#include <vector>
#include <iostream>
#include "ImageSimulator.h"

using namespace cv;
using namespace std;

int main() {
    // Read an image
    string image_path = "imgs/a1.jpg";
    Mat canvas_img = imread(image_path);
    
    int height = canvas_img.size[0];
    int width = canvas_img.size[1];

    // Given a position and a radius, compute a rectangle with a specified image ratio. Then crop parts of the image
    // Allow multiple points & radius (for multiple camera views)
    vector<Point2d> points = { Point2d(800, 800), Point2d(1300, 800), Point2d(1800, 800) };
    vector<float> radius = { 800.0f, 800.0f, 800.0f };

    // Add points and radius into Circle
    vector<Circle*> circles;
    for (int i = 0; i < points.size(); i++) {
        Circle* new_circle = new Circle(i, points[i], radius[i]);
        circles.push_back(new_circle);
    }

    ImageSimulator image_simulator = ImageSimulator(canvas_img, circles);
    image_simulator.load_frames();
    return 0;
}