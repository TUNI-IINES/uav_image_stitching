#include "ImageSimulator.h"
#include "pipeline.h"

ImageSimulator::ImageSimulator(cv::Mat canvas, std::vector<Circle*> circles) :
	canvas_(canvas), circles_(circles){

}

ImageSimulator::~ImageSimulator() {

}

void ImageSimulator::load_frames(int seconds, int fps, int move_x, int move_y, float radius) {
	
	int fourcc = VideoWriter::fourcc('M', 'J', 'P', 'G');
	cv::Size wannabe_size = Size(500, 400);

	VideoWriter outputVideo;
	outputVideo.open("result_video.avi", fourcc, fps, wannabe_size);

	// Get the squares from the circles
	for (int i = 0; i < seconds; i++) {
		// As each second with have a number of frames, call that fps for each image
		for (int f = 0; f < fps; f++) {
			for (Circle* circle : circles_) {
				// Update the circles first
				circle->move(move_x, move_y, radius);
				cv::Mat image = canvas_(cv::Range(circle->y, circle->y + circle->h), cv::Range(circle->x, circle->x + circle->h));
				cv::Mat resized_image;
				cv::resize(image, resized_image, cv::Size(300, 300));
				frames_.push_back(resized_image);
			}
			// Stitch the images here
			std::vector<cv::Mat> src;
			src.push_back(frames_[0]);

			for (int i = 1; i < frames_.size(); i++) {
				src.push_back(frames_[i]);
				Mat stitched_image = pipeline(src);
				src.clear();
				src.push_back(stitched_image);
			}

			cv::resize(src[0], src[0], wannabe_size);

			//cv::imshow("Result stitch", src[0]);
			//cv::waitKey(0);
			outputVideo.write(src[0]);

			// Clear frames
			frames_.clear();
		}
	}

	outputVideo.release();
	cv::destroyAllWindows();
}