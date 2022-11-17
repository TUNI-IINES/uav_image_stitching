#include <vector>
#include "Circle.cpp"

class ImageSimulator
{
public:
	ImageSimulator(cv::Mat canvas);
	~ImageSimulator();

	void init_frames(int fps = 1); 	
	void stitch_frames(std::vector<Circle*> circles, int ms_cnt = 0); 	
	void close_simulator(void); 	

private:
	cv::Mat canvas_;
	std::vector<cv::Mat> frames_;

	cv::Size wannabe_size_;
	cv::VideoWriter outputVideo_;
};

