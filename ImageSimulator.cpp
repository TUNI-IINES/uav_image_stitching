#include "ImageSimulator.h"
#include "pipeline.cpp"

ImageSimulator::ImageSimulator(cv::Mat canvas) :
	canvas_(canvas) {

}

ImageSimulator::~ImageSimulator() {

}

void ImageSimulator::init_frames(int fps) {
	
	int fourcc = VideoWriter::fourcc('M', 'J', 'P', 'G');
	// wannabe_size_ = Size(500, 400);
	wannabe_size_ = Size(1000, 800);
	// outputVideo_.open("result_video.avi", fourcc, fps, wannabe_size_);

}

void ImageSimulator::stitch_frames(std::vector<Circle*> circles, int ms_cnt) {

	for (Circle* circle : circles) {
		// Update the circles first
		cv::Mat image = canvas_(cv::Range(circle->y, circle->y + circle->h), cv::Range(circle->x, circle->x + circle->w));
		cv::Mat resized_image;
		// cv::resize(image, resized_image, cv::Size(300, 300));
		// frames_.push_back(resized_image);
		// Test, if not resizing --> aligning take too long
		// Keep aspect ratio while resizing --> set minimum as 300
		int dim = 300;
		int width = image.cols, height = image.rows;
		width = int(width*dim/height); height = dim;
		//if (width >= height) {
		// 	width = int(width*dim/height); height = dim;
		// } else {
		// 	height = int(height*dim/width); width = dim;
		// }
		cv::resize(image, resized_image, cv::Size(width, height));
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

	std::stringstream sstm;
	sstm << "snapshot_" << ms_cnt << ".jpg";
	imwrite(sstm.str(), src[0]); // Save image before further video resizing
	// cv::imshow("Result stitch", src[0]);
	// cv::waitKey(0);

	// cv::resize(src[0], src[0], wannabe_size_);
	// outputVideo_.write(src[0]);

	// Clear frames
	frames_.clear();
}

void ImageSimulator::close_simulator(void) {
	// outputVideo_.release();
	cv::destroyAllWindows();
}