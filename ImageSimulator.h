#include <vector>
#include "Circle.h"

class ImageSimulator
{
public:
	ImageSimulator(cv::Mat canvas, std::vector<Circle*> circles);
	~ImageSimulator();

	void load_frames(int seconds = 30, int fps = 1, int move_x = 10, int move_y = 10, float radius = -1.0f);

private:
	cv::Mat canvas_;
	std::vector<Circle*> circles_;
	std::vector<cv::Mat> frames_;
};

