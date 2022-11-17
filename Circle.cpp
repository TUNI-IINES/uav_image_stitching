#include "Circle.h"

Circle::Circle(int id, cv::Point2d coords, float radius) :
	id_(id),
	coords_(coords),
	radius_(radius) {
	x = 0;
	y = 0;
	h = 0;
	update_square();
}

Circle::~Circle() {

}

void Circle::move(int move_x, int move_y, float radius=-1) {
	
	if (coords_.x + move_x - radius_ >= 0) {
		coords_.x += move_x;
	}
	if (coords_.y + move_y - radius_ >= 0) {
		coords_.y += move_y;
	}
	if (radius != -1) {
		radius_ = radius;
	}
	update_square();
}

void Circle::update_square() {
	// Calculate the square to crop from a big canvas
	if (coords_.x - radius_ >= 0) {
		x = coords_.x - radius_;
	}
	else {
		x = 0;
	}
	if (coords_.y - radius_ >= 0) {
		y = coords_.y - radius_;
	}
	else {
		y = 0;
	}
	h = radius_ * 2;
}