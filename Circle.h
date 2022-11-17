#include "opencv2/opencv.hpp"

class Circle
{
public:
	Circle(int id, cv::Point2d coords, float radius);
	~Circle();

	void move(int move_x, int move_y, float radius);
	void update_pose(int pos_x, int pos_y, float radius);
	int x, y, w, h;

private:
	void update_square();
	int id_;
	// Circle center
	cv::Point2d coords_;
	float radius_x;
	float radius_y;
};

