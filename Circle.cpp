#include "Circle.h"

Circle::Circle(int id, cv::Point2d coords, float radius) :
	id_(id),
	coords_(coords),
	radius_y(radius) {
	x = 0;
	y = 0;
	w = 0;
	h = 0;
	update_square();
}

Circle::~Circle() {

}

void Circle::move(int move_x, int move_y, float radius=-1) {
	
	if (coords_.x + move_x - radius_x >= 0) {
		coords_.x += move_x;
	}
	if (coords_.y + move_y - radius_y >= 0) {
		coords_.y += move_y;
	}
	if (radius != -1) {
		radius_y = radius;
		radius_x = radius*4/3;
	}
	update_square();
}

void Circle::update_pose(int pos_x, int pos_y, float radius=-1) {
	
	if (pos_x - radius_x >= 0) {
		coords_.x = pos_x;
	}
	if (pos_y - radius_y >= 0) {
		coords_.y = pos_y;
	}
	if (radius != -1) {
		radius_y = radius;
		radius_x = radius*4/3;
	}
	update_square();
}


void Circle::update_square() {
	// Calculate the square to crop from a big canvas
	if (coords_.x - radius_x >= 0) {
		x = coords_.x - radius_x;
	}
	else {
		x = 0;
	}
	if (coords_.y - radius_y >= 0) {
		y = coords_.y - radius_y;
	}
	else {
		y = 0;
	}
	h = radius_y * 2;
	w = radius_x * 2;
}