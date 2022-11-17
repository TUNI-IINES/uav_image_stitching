#pragma once
#ifndef M_WATERSHED_H__
#define M_WATERSHED_H__
//#include "stdafx.h"
#include <vector>
#include <queue>
#include <algorithm>
#include <time.h>
#include <math.h>
#include <iostream>
#include <set>
#include <map>
#include <string.h>
#include <float.h>
#include <stdlib.h>
#include <stack>
#include <fstream>

using namespace std;

const int INIT = -3;
const int WATERSHED = -1;
const int MAXN = 8000;
const int lmax = 0x7FFFFFFF;

#define Num_of_Queues  750
#define inf FLT_MAX

typedef struct {
	int x;
	int y;
}zwPoint;

typedef struct iPoint {
	int x;
	int y;
}spoint;

typedef struct SuperPixell {
	float RegionChangeRate;
	float Avg;
	int PixelNum;
	spoint* pixel;
	float SuperPixelDifferece;
}SuperPixel;

double factorial(double num);
void create_filters(int WinSize, double* smofil_, double* diffil_);
void gaussian_diff_filter(double* pImage, float* pGrad, int ImgSizeX, int ImgSizeY, int WinSize);
void get_gaussian_gradient(double* R, double* G, double* B, int iWidth, int iHeight, int WinSize, float* deltar);
void set_uniform_seeds(int* Label, int* mask, float* deltar, int iWidth, int iHeight, int WinSize, std::vector< zwPoint >& seedPos);
void ExpandBasin(int center, int pos, int* Label, int* area);
void superpixel_flood(double lambda, std::vector< zwPoint > m_SeedPos, float* deltar, int WinSize, int ImgSizeX, int ImgSizeY, int* Label, int* mask);
void findedges(int* Label, unsigned char* edge, int iWidth, int iHeight);
double distance(double* R, double* G, double* B, int pos1, int pos2);
void ProcEdgePixs(int* Label, int* mask, double* R, double* G, double* B, std::vector< zwPoint > seedPos, int ImgSizeX, int ImgSizeY);
void SmootBoundaries(int* Label, int ImgSizeX, int ImgSizeY);
void scw_rgb(double* R, double* G, double* B, int* Label, int* mask, int iWidth, int iHeight, int WinSize, double lambda);
void classify(int* Label, int* classmask, int ImgsizeX1, int ImgsizeY1, int ImgsizeX2, int ImgsizeY2, int ImgsizeX, int ImgsizeY, int& num1, int& num2);
int cmp(const void* a, const void* b);
void BuildRegionAdjacencyMatrix(int** RegionAdjacencyMatrix, int* Label, int row, int col);
void GetGradient(double* picture1, double* picture2, float** Gradient, int sizeX, int sizeY, int SuperpixelNum, SuperPixel** supPixel);
void GetSuperpixelFeature(double* picture1, double* picture2, int* Label, int sizeX, int sizeY, int SuperpixelNum, SuperPixel** supPixel);
void release_superpixel(int SuperpixelNum, SuperPixel** supPixel);
void BuildCostWeightMatrix(int** RegionAdjacencyMatrix, float** CostWeightMatrix, int MaxLabel, SuperPixel** supPixel, float weight1, float weight2, float weight3);
void init(queue<int>& q, int* Path);
float Min(float a, float b);
float bfs(queue<int>& q, float** Net, int* Path, float* Lv, int& superpixel_num, int& Start, int& End);
void cutpart(float** Net, int* Path, int& pre_superpixelNum, int* out, int& t);
void Ford_Fulkerson(queue<int>& q, float** Net, float** Net1, int* Path, float* Lv, int& superpixel_num, int& nn, int& Start, int& End, int& count1, int& count2, int* c1, int* c2, int& pre_superpixelNum, int* out, int& t, float** data);
void mincut(int superpixel_num, int count1, int count2, float** data, int* c1, int* c2, int* out);
void dfs_class(int* cutline, int* Label, int* symbol, int* Labelclass, int curpos, int iWidth, int iHeight);
std::pair<int*, int*> seamcut(int iWidth, int iHeight, int* mask, int* Label, int MaxLabel, std::vector<float> cost, int ImgsizeX1, int ImgsizeY1, int ImgsizeX2, int ImgsizeY2);


class color_node {
private:
	double* sColor;
	int m_len;
public:
	color_node() { sColor = new double[m_len]; }
	~color_node() { delete[] sColor; }
};

#endif 
