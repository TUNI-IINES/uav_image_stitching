// seamcut.cpp : 
/**********************mincut****************************************/
//#include "stdafx.h"
#include"seamcut.h"
#include "stdint.h"


// function	factorial
double factorial(double num)
{
	if (num == 0 || num == 1)
		return 1;
	return (num * factorial(num - 1));
}

// function	create_filters
void create_filters(int WinSize, double* smofil_, double* diffil_)
{
	int i;
	double w;
	int WL_ = WinSize / 2;
	for (i = -WL_; i <= WL_; i++)
	{
		w = exp(double(-2 * WL_)) * factorial(2 * WL_) / (factorial(WL_ - i) * factorial(WL_ + i));
		smofil_[i + WL_] = w;
		diffil_[i + WL_] = (2 * i * w) / WL_;
	}
	double t1 = 0.0f;
	double t2 = 0.0f;
	for (i = 0; i < WinSize; i++)
	{
		t1 = smofil_[i] > t1 ? smofil_[i] : t1;
		t2 = diffil_[i] > t2 ? diffil_[i] : t2;
	}
	for (i = 0; i < WinSize; i++)
	{
		smofil_[i] /= t1;
		diffil_[i] /= t2;
	}
}

// function	gaussian_diff_filter
void gaussian_diff_filter(double* pImage, float* pGrad, int ImgSizeX, int ImgSizeY, int WinSize)
{

	//smooth filter
	double* sf = new double[WinSize];
	//diff filter
	double* df = new double[WinSize];
	//unsigned char* im;
	double* tim;
	double sum = 0;
	int i, j, k;
	//create kernels
	create_filters(WinSize, sf, df);
	int m_nWidth = ImgSizeX;
	int m_nHeight = ImgSizeY;
	float* pGradX = new float[m_nWidth * m_nHeight];
	float* pGradY = new float[m_nWidth * m_nHeight];
	//im = cim->im_;
	tim = new double[m_nWidth * m_nHeight];
	for (i = 0; i < m_nWidth * m_nHeight; i++)
	{
		pGradX[i] = pGradY[i] = 0;
		tim[i] = pImage[i];
	}
	int WL_ = WinSize / 2;
	//filter image x
	//smooth on y
	for (i = 0; i < m_nWidth; i++)
	{
		for (j = WL_; j < (m_nHeight - WL_); j++)
		{
			sum = 0;
			for (k = -WL_; k <= WL_; k++)
				sum += (float)sf[k + WL_] * pImage[(j + k) * m_nWidth + i];
			tim[j * m_nWidth + i] = sum;
		}
	}
	//diff on x
	for (j = 0; j < m_nHeight; j++)
	{
		for (i = WL_; i < (m_nWidth - WL_); i++)
		{
			sum = 0;
			for (k = -WL_; k <= WL_; k++)
				sum += df[k + WL_] * tim[j * m_nWidth + i + k];
			pGradX[j * m_nWidth + i] = (float)(sum);
		}
	}
	//filter image y
	for (i = 0; i < m_nWidth * m_nHeight; i++)
		tim[i] = pImage[i];
	//  im = cim->im_;
	//smooth on x
	for (j = 0; j < m_nHeight; j++)
	{
		for (i = WL_; i < (m_nWidth - WL_); i++)
		{
			sum = 0;
			for (k = -WL_; k <= WL_; k++)
				sum += (float)sf[k + WL_] * pImage[j * m_nWidth + i + k];
			tim[j * m_nWidth + i] = sum;
		}
	}
	//diff on y
	for (i = 0; i < m_nWidth; i++)
	{
		for (j = WL_; j < (m_nHeight - WL_); j++)
		{
			sum = 0;
			for (k = -WL_; k <= WL_; k++)
				sum += df[k + WL_] * tim[(j + k) * m_nWidth + i];
			pGradY[j * m_nWidth + i] = (float)(sum);
		}
	}
	double weight = 1.0 / (0.3333333 * WinSize);
	double tmp;
	for (i = 0; i < m_nWidth * m_nHeight; i++)
	{
		tmp = pow((double)(pGradX[i] * pGradX[i] + pGradY[i] * pGradY[i]), 0.5);
		pGrad[i] = (float)(weight * tmp < 255 ? weight * tmp : 255);
	}
	delete[]pGradX;
	delete[]pGradY;
	delete[]tim;
	delete[]sf;
	delete[]df;
}

// function	get_gaussian_gradient
void get_gaussian_gradient(double* R, double* G, double* B, int iWidth, int iHeight, int WinSize, float* deltar)
{
	float* deltar_r = new float[iWidth * iHeight];
	float* deltar_g = new float[iWidth * iHeight];
	float* deltar_b = new float[iWidth * iHeight];
	gaussian_diff_filter(R, deltar_r, iWidth, iHeight, WinSize);
	gaussian_diff_filter(G, deltar_g, iWidth, iHeight, WinSize);
	gaussian_diff_filter(B, deltar_b, iWidth, iHeight, WinSize);
	int i;
	for (i = 0; i < iWidth * iHeight; i++)
	{
		deltar[i] = (deltar_r[i] + deltar_g[i] + deltar_b[i]) * 0.33333f * 255;
	}
	delete[]deltar_r; deltar_r = NULL;
	delete[]deltar_g; deltar_g = NULL;
	delete[]deltar_b; deltar_b = NULL;
}

// function	set_uniform_seeds
void set_uniform_seeds(int* Label, int* mask, float* deltar, int iWidth, int iHeight, int WinSize, std::vector< zwPoint >& seedPos)
{
	int i;
	int tmppos;
	int seedpos;
	int _WinSize = WinSize;
	int imagelen = iWidth * iHeight;
	if (WinSize < 4)
	{
		WinSize = 4;
	}
	for (i = 0; i < imagelen; i++)
	{
		Label[i] = INIT;
	}
	//Set seed points evenly in the image
	int x, y;
	float minGradient = 1000;
	int blockSizeX = iWidth / _WinSize;
	float xDist = ((float)iWidth) / blockSizeX;
	int blockSizeY = iHeight / _WinSize;
	float yDist = ((float)iHeight) / blockSizeY;
	int j;
	int CentroX, CentroY;
	zwPoint tmpFPoint;
	int SuperNum = 0;
	for (j = 0; j < blockSizeY; j++)
	{
		for (i = 0; i < blockSizeX; i++)
		{
			CentroY = (int)((0.5 + j) * yDist);
			CentroX = (int)((0.5 + i) * xDist);
			minGradient = 1000;
			if (mask[CentroY * iWidth + CentroX] == 1)
			{
				for (y = CentroY - 2; y < CentroY + 2; y++)
				{
					for (x = CentroX - 2; x < CentroX + 2; x++)
					{
						tmppos = y * iWidth + x;
						if (deltar[tmppos] < minGradient)
						{
							minGradient = deltar[tmppos];
							seedpos = tmppos;
						}
					}
				}
				Label[CentroY * iWidth + CentroX] = SuperNum;
				tmpFPoint.y = seedpos / iWidth;
				tmpFPoint.x = seedpos - iWidth * tmpFPoint.y;
				//Record the location of the seed point
				seedPos.push_back(tmpFPoint);
				SuperNum = SuperNum + 1;
			}
		}
	}
}

// function	ExpandBasin
void ExpandBasin(int center, int pos, int* Label, int* area)
{
	if (Label[pos] > WATERSHED)
	{
		//Neither a watershed nor an area
		if (Label[center] < WATERSHED)
		{
			//Set it as the area of the neighboring point
			Label[center] = Label[pos];
			area[Label[pos]] += 1;
		}
		else if (Label[center] == WATERSHED)
		{

		}
		//Belongs to a certain area or is inconsistent with the adjacent pixel area label
		else if (Label[center] != Label[pos])
		{
			//The original district is different from the current district, which is set as a watershed
			Label[center] = WATERSHED;
		}
	}
}

// function	superpixel_flood
void superpixel_flood(double lambda, std::vector< zwPoint > m_SeedPos, float* deltar, int WinSize, int ImgSizeX, int ImgSizeY, int* Label, int* mask)
{


	int i;
	int Imglen = ImgSizeX * ImgSizeY;
	//the first and last row
	for (i = 0; i < ImgSizeX; i++)
	{
		Label[i] = -1;
		Label[Imglen - 1 - i] = -1;
	}
	//the first and last column
	for (i = 0; i < ImgSizeY; i++)
	{
		Label[i * ImgSizeX] = -1;
		Label[(i + 1) * ImgSizeX - 1] = -1;
	}
	double  mean = 0;
	for (i = 0; i < Imglen; i++)
	{
		mean += deltar[i];
	}
	mean = mean / Imglen;
	double beta;
	beta = pow(mean, (double)0.8f);
	double m_bShapeConstraint[Num_of_Queues];
	memset(m_bShapeConstraint, 0, Num_of_Queues * sizeof(double));
	int ObjNum = (int)m_SeedPos.size();
	int* area = new int[ObjNum];
	memset(area, 0, sizeof(int) * ObjNum);
	for (i = 0; i < Num_of_Queues; i++)
	{
		m_bShapeConstraint[i] = 1.0f / (1.0f + pow((double)(i / beta), (int)2)) * lambda;
	}
	double* distweight = new double[500];
	for (i = 0; i < 500; i++)
	{
		double tmp = -i;
		distweight[i] = 1 - exp(tmp * 0.2);
	}
	unsigned char* isProcessed = new unsigned char[ImgSizeX * ImgSizeY];
	memset(isProcessed, 0, ImgSizeX * ImgSizeY * sizeof(unsigned char));
	std::queue< zwPoint > PriQue[Num_of_Queues];
	int j;
	int TmpY, TmpPos;
	zwPoint TmpPt;
	for (i = 1; i < ImgSizeY - 1; i++)
	{
		TmpY = i * ImgSizeX;
		for (j = 1; j < ImgSizeX - 1; j++)
		{
			TmpPos = TmpY + j;
			if (Label[TmpPos] < WATERSHED)
				if (Label[TmpPos - 1] > WATERSHED
					|| Label[TmpPos + 1] > WATERSHED
					|| Label[TmpPos + ImgSizeX] > WATERSHED
					|| Label[TmpPos - ImgSizeX] > WATERSHED
					)
				{
					TmpPt.x = j; TmpPt.y = i;
					PriQue[0].push(TmpPt);
					isProcessed[TmpPos] = 1;
				}
		}
	}
	zwPoint TmpPt2;
	unsigned short Index[Num_of_Queues];
	for (i = 0; i < Num_of_Queues; i++)
	{
		Index[i] = i;
	}
	for (i = 0; i < Num_of_Queues; i++)
	{

		int k = 0;
		for (k = 0; k < i; k++)
		{
			Index[k] = i;
		}
		double dis2seed = 0.0f;
		zwPoint seedPos;
		while (!PriQue[i].empty())
		{
			TmpPt = PriQue[i].front();
			PriQue[i].pop();
			TmpPos = TmpPt.y * ImgSizeX + TmpPt.x;
			ExpandBasin(TmpPos, TmpPos - 1, Label, area);        //left
			ExpandBasin(TmpPos, TmpPos + 1, Label, area);        //right
			ExpandBasin(TmpPos, TmpPos - ImgSizeX, Label, area); //up
			ExpandBasin(TmpPos, TmpPos + ImgSizeX, Label, area); //down
			if (Label[TmpPos] >= 0)
			{
				seedPos.y = m_SeedPos[Label[TmpPos]].y;
				seedPos.x = m_SeedPos[Label[TmpPos]].x;
				if (Label[TmpPos - 1] < WATERSHED && !isProcessed[TmpPos - 1])
				{
					TmpPt2.x = TmpPt.x - 1; TmpPt2.y = TmpPt.y;
					dis2seed = (seedPos.x - TmpPt2.x) * (seedPos.x - TmpPt2.x) + (seedPos.y - TmpPt2.y) * (seedPos.y - TmpPt2.y);
					dis2seed = sqrt((double)dis2seed);
					double p = Index[int(deltar[TmpPos - 1] * distweight[int(dis2seed)])] + m_bShapeConstraint[i] * dis2seed;
					PriQue[(int)p].push(TmpPt2);
					isProcessed[TmpPos - 1] = 1;
				}
				if (Label[TmpPos + 1] < WATERSHED && !isProcessed[TmpPos + 1])
				{
					TmpPt2.x = TmpPt.x + 1; TmpPt2.y = TmpPt.y;
					dis2seed = (seedPos.x - TmpPt2.x) * (seedPos.x - TmpPt2.x) + (seedPos.y - TmpPt2.y) * (seedPos.y - TmpPt2.y);
					dis2seed = sqrt((double)dis2seed);
					double p = Index[int(deltar[TmpPos + 1] * distweight[int(dis2seed)])] + m_bShapeConstraint[i] * dis2seed;
					PriQue[(int)p].push(TmpPt2);
					isProcessed[TmpPos + 1] = 1;
				}
				if (Label[TmpPos - ImgSizeX] < WATERSHED && !isProcessed[TmpPos - ImgSizeX])
				{
					TmpPt2.x = TmpPt.x; TmpPt2.y = TmpPt.y - 1;
					dis2seed = (seedPos.x - TmpPt2.x) * (seedPos.x - TmpPt2.x) + (seedPos.y - TmpPt2.y) * (seedPos.y - TmpPt2.y);
					dis2seed = sqrt((double)dis2seed);
					double p = Index[int(deltar[TmpPos - ImgSizeX] * distweight[int(dis2seed)])] + m_bShapeConstraint[i] * dis2seed;
					PriQue[(int)p].push(TmpPt2);
					isProcessed[TmpPos - ImgSizeX] = 1;
				}
				if (Label[TmpPos + ImgSizeX] < WATERSHED && !isProcessed[TmpPos + ImgSizeX])
				{
					TmpPt2.x = TmpPt.x; TmpPt2.y = TmpPt.y + 1;
					dis2seed = (seedPos.x - TmpPt2.x) * (seedPos.x - TmpPt2.x) + (seedPos.y - TmpPt2.y) * (seedPos.y - TmpPt2.y);
					dis2seed = sqrt((double)dis2seed);
					double p = Index[int(deltar[TmpPos + ImgSizeX] * distweight[int(dis2seed)])] + m_bShapeConstraint[i] * dis2seed;
					PriQue[(int)p].push(TmpPt2);
					isProcessed[TmpPos + ImgSizeX] = 1;
				}
			}
		}
	}
	delete[]area;
	area = NULL;
	delete[]isProcessed;
	isProcessed = NULL;
	delete[]distweight;
	distweight = NULL;
}

// function	findedges
void findedges(int* Label, unsigned char* edge, int iWidth, int iHeight)
{
	int i;
	for (i = 0; i < iHeight * iWidth; i++)
	{
		if (Label[i] < 0)
		{
			edge[i] = 255;
		}
		else edge[i] = 0;
	}
}

// function	distance
double distance(double* R, double* G, double* B, int pos1, int pos2)
{
	double dist = 0.0f;
	dist += (R[pos1] - R[pos2]) * (R[pos1] - R[pos2]);
	dist += (G[pos1] - R[pos2]) * (G[pos1] - R[pos2]);
	dist += (B[pos1] - R[pos2]) * (B[pos1] - R[pos2]);
	return dist;
}

// function	ProcEdgePixs
void ProcEdgePixs(int* Label, int* mask, double* R, double* G, double* B, std::vector< zwPoint > seedPos, int ImgSizeX, int ImgSizeY)
{
	int x[4] = { -1,1,-ImgSizeX,ImgSizeX };
	int i, j;
	int xstart, tmppos;
	unsigned char* edge = new unsigned char[ImgSizeX * ImgSizeY];
	findedges(Label, edge, ImgSizeX, ImgSizeY);
	int k;
	double Dist;
	double min = 10000;
	int local;
	for (i = 1; i < ImgSizeY - 1; i++)
	{
		xstart = i * ImgSizeX;
		for (j = 1; j < ImgSizeX - 1; j++)
		{
			tmppos = xstart + j;
			if (edge[tmppos] == 255)
			{
				min = 999999;
				local = 0;
				for (k = 0; k < 4; k++)
				{
					if (edge[tmppos + x[k]] != 255)
					{
						Dist = distance(R, G, B, tmppos, tmppos + x[k]);
						Dist += sqrt((double)((seedPos[Label[tmppos + x[k]]].x - j) * (seedPos[Label[tmppos + x[k]]].x - j) + (seedPos[Label[tmppos + x[k]]].y - i) * (seedPos[Label[tmppos + x[k]]].y - i)));
						if (Dist < min)
						{
							local = k;
							min = Dist;
						}
					}
				}
				Label[tmppos] = Label[tmppos + x[local]];
				if (Label[tmppos] < 0)
				{
					for (k = 0; k < 4; k++)
					{
						if (Label[tmppos + x[k]] >= 0)
						{
							Label[tmppos] = Label[tmppos + x[k]];
							break;
						}
					}
				}
			}
		}
	}
	for (i = 0; i < ImgSizeY * ImgSizeX; i++)
	{
		if (mask[i] == 0)
			Label[i] = -1;
	}
	delete[]edge; edge = NULL;
}

// function	SmootBoundaries
void SmootBoundaries(int* Label, int ImgSizeX, int ImgSizeY)
{
	int up, down, left, right;
	int flag = 0;
	int m_b[4];
	int i, j, xstart, tmppos;
	for (i = 1; i < ImgSizeY - 1; i++)
	{
		xstart = i * ImgSizeX;
		for (j = 1; j < ImgSizeX - 1; j++)
		{
			tmppos = xstart + j;
			up = tmppos - ImgSizeX;
			down = tmppos + ImgSizeX;
			left = tmppos - 1;
			right = tmppos + 1;
			flag = 0;
			memset(m_b, 0, 4 * sizeof(int));
			if (Label[tmppos] == Label[up]) { flag++; m_b[0] = 1; }
			if (Label[tmppos] == Label[left]) { flag++; m_b[1] = 1; }
			if (Label[tmppos] == Label[down]) { flag++; m_b[2] = 1; }
			if (Label[tmppos] == Label[right]) { flag++; m_b[3] = 1; }
			if (flag == 1)
			{
				if (Label[up] == Label[down] || Label[left] == Label[right])
				{
					if (!m_b[0]) { Label[tmppos] = Label[up]; }
					else if (!m_b[1]) { Label[tmppos] = Label[left]; }
					else if (!m_b[2]) { Label[tmppos] = Label[down]; }
					else if (!m_b[3]) { Label[tmppos] = Label[right]; }
				}
			}
		}
	}
}

// function	scw_rgb
void scw_rgb(double* R, double* G, double* B, int* Label, int* mask, int iWidth, int iHeight, int WinSize, double lambda)
{

	float* deltar = new float[iWidth * iHeight];
	get_gaussian_gradient(R, G, B, iWidth, iHeight, 3, deltar);
	std::vector< zwPoint > seedPos;
	set_uniform_seeds(Label, mask, deltar, iWidth, iHeight, WinSize, seedPos);
	superpixel_flood(lambda, seedPos, deltar, WinSize, iWidth, iHeight, Label, mask);

	delete[]deltar; deltar = NULL;
	ProcEdgePixs(Label, mask, R, G, B, seedPos, iWidth, iHeight);
	SmootBoundaries(Label, iWidth, iHeight);
}

// function	classify
void classify(int* Label, int* classmask, int ImgsizeX1, int ImgsizeY1, int ImgsizeX2, int ImgsizeY2, int ImgsizeX, int ImgsizeY, int& num1, int& num2)
{
	int* flag = new int[ImgsizeX * ImgsizeY];
	memset(flag, 0, sizeof(int) * ImgsizeX * ImgsizeY);
	int* margin = new int[ImgsizeX * ImgsizeY];
	memset(margin, 0, sizeof(int) * ImgsizeX * ImgsizeY);
	int large, little;
	if ((ImgsizeY1 * ImgsizeX + ImgsizeX1) < (ImgsizeY2 * ImgsizeX + ImgsizeX2))
	{
		large = ImgsizeY2 * ImgsizeX + ImgsizeX2;
		little = ImgsizeY1 * ImgsizeX + ImgsizeX1;
	}
	else
	{
		large = ImgsizeY1 * ImgsizeX + ImgsizeX1;
		little = ImgsizeY2 * ImgsizeX + ImgsizeX2;
	}
	int direction1[4] = { -1,-ImgsizeX,1,ImgsizeX };
	for (int i = 0; i < ImgsizeX * ImgsizeY; i++)
	{
		for (int k = 0; k < 4; k++)
			if (Label[i] >= 0 && i + direction1[k] >= 0 && i + direction1[k] < ImgsizeX * ImgsizeY && Label[i + direction1[k]] < 0)
				margin[i] = 1;
	}
	int direction2[8] = { -ImgsizeX - 1,-1,ImgsizeX - 1,ImgsizeX,ImgsizeX + 1,1,-ImgsizeX + 1,-ImgsizeX };
	int curPos = little;
	bool reach = false;
	//class1
	flag[little] = 1;
	stack<int> t1;
	while (!t1.empty())
		t1.pop();
	t1.push(curPos);
	while (!t1.empty())
	{
		int temp1 = t1.top();
		t1.pop();
		if (temp1 == large)
			break;
		for (int k = 0; k < 8; k++)
		{
			if (temp1 + direction2[k] >= 0 && temp1 + direction2[k] < ImgsizeX * ImgsizeY && temp1 + direction2[k] == large)
				reach = true;
		}
		if (reach)
			continue;
		else
			for (int k = 0; k < 8; k++)
			{
				if (temp1 + direction2[k] >= 0 && temp1 + direction2[k] < ImgsizeX * ImgsizeY && margin[temp1 + direction2[k]] == 1 && flag[temp1 + direction2[k]] == 0)
				{
					flag[temp1 + direction2[k]] = 1;
					classmask[temp1 + direction2[k]] = 1;
					num1++;
					t1.push(temp1 + direction2[k]);
				}
			}
	}
	flag[little] = 0;
	flag[large] = 1;
	reach = false;
	//class2
	curPos = large;
	stack<int> t2;
	while (!t2.empty())
		t2.pop();
	t2.push(curPos);
	while (!t2.empty())
	{
		int temp2 = t2.top();
		t2.pop();
		if (temp2 == little)
			break;
		for (int k = 0; k < 8; k++)
		{
			if (temp2 + direction2[k] >= 0 && temp2 + direction2[k] < ImgsizeX * ImgsizeY && margin[temp2 + direction2[k]] == 1 && flag[temp2 + direction2[k]] == 0)
			{
				flag[temp2 + direction2[k]] = 1;
				classmask[temp2 + direction2[k]] = 2;
				num2++;
				t2.push(temp2 + direction2[k]);
			}
		}
	}
	for (int k = 0; k < 8; k++)
	{
		if (little + direction2[k] >= 0 && little + direction2[k] < ImgsizeX * ImgsizeY && little + 2 * direction2[k] >= 0 && little + 2 * direction2[k] < ImgsizeX * ImgsizeY && margin[little + direction2[k]] == 1 && flag[little + direction2[k]] == 1 && margin[little + 2 * direction2[k]] == 1 && flag[little + 2 * direction2[k]] == 1)
			classmask[little + direction2[k]] = classmask[little + 2 * direction2[k]];
	};
	delete[]flag; flag = NULL;
	delete[]margin; margin = NULL;
}

// function	cmp
int cmp(const void* a, const void* b)
{
	if (*(float*)a >= *(float*)b) return 1;
	else return -1;
}

// function	BuildRegionAdjacencyMatrix
void BuildRegionAdjacencyMatrix(int** RegionAdjacencyMatrix, int* Label, int row, int col)
{
	for (int i = 1; i < row - 1; i++)
	{
		for (int j = 1; j < col - 1; j++)
		{
			if (Label[i * col + j] >= 0)
			{
				if ((Label[(i - 1) * col + j - 1] >= 0) && (Label[i * col + j] != Label[(i - 1) * col + j - 1]))
				{
					RegionAdjacencyMatrix[Label[i * col + j]][Label[(i - 1) * col + j - 1]] = 1;
					RegionAdjacencyMatrix[Label[(i - 1) * col + j - 1]][Label[i * col + j]] = 1;
				}
				if ((Label[(i - 1) * col + j] >= 0) && (Label[i * col + j] != Label[(i - 1) * col + j]))
				{
					RegionAdjacencyMatrix[Label[i * col + j]][Label[(i - 1) * col + j]] = 1;
					RegionAdjacencyMatrix[Label[(i - 1) * col + j]][Label[i * col + j]] = 1;
				}
				if ((Label[(i - 1) * col + j + 1] >= 0) && (Label[i * col + j] != Label[(i - 1) * col + j + 1]))
				{
					RegionAdjacencyMatrix[Label[i * col + j]][Label[(i - 1) * col + j + 1]] = 1;
					RegionAdjacencyMatrix[Label[(i - 1) * col + j + 1]][Label[i * col + j]] = 1;
				}
				if ((Label[i * col + j - 1] >= 0) && (Label[i * col + j] != Label[i * col + j - 1]))
				{
					RegionAdjacencyMatrix[Label[i * col + j]][Label[i * col + j - 1]] = 1;
					RegionAdjacencyMatrix[Label[i * col + j - 1]][Label[i * col + j]] = 1;
				}
				if ((Label[i * col + j + 1] >= 0) && (Label[i * col + j] != Label[i * col + j + 1]))
				{
					RegionAdjacencyMatrix[Label[i * col + j]][Label[i * col + j + 1]] = 1;
					RegionAdjacencyMatrix[Label[i * col + j + 1]][Label[i * col + j]] = 1;
				}
				if ((Label[(i + 1) * col + j - 1] >= 0) && (Label[i * col + j] != Label[(i + 1) * col + j - 1]))
				{
					RegionAdjacencyMatrix[Label[i * col + j]][Label[(i + 1) * col + j - 1]] = 1;
					RegionAdjacencyMatrix[Label[(i + 1) * col + j - 1]][Label[i * col + j]] = 1;
				}
				if ((Label[(i + 1) * col + j] >= 0) && (Label[i * col + j] != Label[(i + 1) * col + j]))
				{
					RegionAdjacencyMatrix[Label[i * col + j]][Label[(i + 1) * col + j]] = 1;
					RegionAdjacencyMatrix[Label[(i + 1) * col + j]][Label[i * col + j]] = 1;
				}
				if ((Label[(i + 1) * col + j + 1] >= 0) && (Label[i * col + j] != Label[(i + 1) * col + j + 1]))
				{
					RegionAdjacencyMatrix[Label[i * col + j]][Label[(i + 1) * col + j + 1]] = 1;
					RegionAdjacencyMatrix[Label[(i + 1) * col + j + 1]][Label[i * col + j]] = 1;
				}
			}
		}
	}
}

// function	GetGradient
void GetGradient(double* picture1, double* picture2, float** Gradient, int sizeX, int sizeY, int SuperpixelNum, SuperPixel** supPixel)
{
	int i, j;
	float** subpicture = (float**)malloc(sizeof(float*) * sizeY);
	for (i = 0; i < sizeY; i++)
		subpicture[i] = (float*)malloc(sizeof(float) * sizeX);
	float** GradientX = (float**)malloc(sizeof(float*) * sizeY);
	for (i = 0; i < sizeY; i++)
		GradientX[i] = (float*)malloc(sizeof(float) * sizeX);
	float** GradientY = (float**)malloc(sizeof(float*) * sizeY);
	for (i = 0; i < sizeY; i++)
		GradientY[i] = (float*)malloc(sizeof(float) * sizeX);
	for (i = 0; i < sizeY; i++)
	{
		for (j = 0; j < sizeX; j++)
			subpicture[i][j] = (float)(picture1[i * sizeX + j] - picture2[i * sizeX + j]);
	}
	for (i = 0; i < sizeY; i++)
	{
		GradientX[i][0] = 0;
		for (j = 1; j < sizeX; j++)
			GradientX[i][j] = subpicture[i][j] - subpicture[i][j - 1];
	}
	for (i = 0; i < sizeX; i++)
	{
		GradientY[0][i] = 0;
		for (j = 1; j < sizeY; j++)
			GradientY[j][i] = subpicture[j][i] - subpicture[j - 1][i];
	}
	for (i = 0; i < sizeY; i++)
	{
		for (j = 1; j < sizeX; j++)
			Gradient[i][j] = sqrt(GradientX[i][j] * GradientX[i][j] + GradientY[i][j] * GradientY[i][j]);
	}
	for (i = 0; i < sizeY; i++)
	{
		free(GradientX[i]);
		free(GradientY[i]);
		free(subpicture[i]);
	}
	free(GradientX);
	free(GradientY);
	free(subpicture);
}

// function	GetSuperpixelFeature
void GetSuperpixelFeature(double* picture1, double* picture2, int* Label, int sizeX, int sizeY, int SuperpixelNum, SuperPixel** supPixel)
{
	int i, j, k, l;
	int thisPixelNum = 0;
	for (k = 0; k < SuperpixelNum; k++)
	{
		for (i = 0; i < sizeY; i++)
		{
			for (j = 0; j < sizeX; j++)
				if (Label[i * sizeX + j] == k)   thisPixelNum++;
		}
		supPixel[k]->pixel = (spoint*)malloc(sizeof(spoint) * thisPixelNum);
		thisPixelNum = 0;
		for (i = 0; i < sizeY; i++)
		{
			for (j = 0; j < sizeX; j++)
				if (Label[i * sizeX + j] == k)
				{
					spoint* tempspspoint;
					tempspspoint = supPixel[k]->pixel;
					spoint* tempPoint = &tempspspoint[thisPixelNum];
					tempPoint->x = i;
					tempPoint->y = j;
					thisPixelNum++;
				}
		}
		float sum = 0;
		supPixel[k]->PixelNum = thisPixelNum;
		int tPixelNum = thisPixelNum;
		for (thisPixelNum = 0; thisPixelNum < tPixelNum; thisPixelNum++)
		{
			int pictureX, pictureY;
			spoint* tempspspoint;
			tempspspoint = supPixel[k]->pixel;
			spoint tempPoint = tempspspoint[thisPixelNum];
			pictureX = tempPoint.x;
			pictureY = tempPoint.y;
			sum += (float)picture1[pictureX * sizeX + pictureY];
		}
		supPixel[k]->Avg = sum / tPixelNum;
	}
	float** Ncc = (float**)malloc(sizeof(float*) * sizeY);
	for (i = 0; i < sizeY; i++)
		Ncc[i] = (float*)malloc(sizeof(float) * sizeX);
	float* PD = (float*)malloc(sizeof(float) * sizeX * sizeY);
	float** newPicture1 = (float**)malloc(sizeof(float*) * (4 + sizeY));
	for (i = 0; i < sizeY + 4; i++)
		newPicture1[i] = (float*)malloc(sizeof(float) * (sizeX + 4));
	for (i = 0; i < sizeX + 4; i++)
		newPicture1[0][i] = newPicture1[1][i] = newPicture1[sizeY + 2][i] = newPicture1[sizeY + 3][i] = -1;
	for (i = 0; i < sizeY + 4; i++)
		newPicture1[i][0] = newPicture1[i][1] = newPicture1[i][sizeX + 2] = newPicture1[i][sizeX + 3] = -1;
	for (i = 0; i < sizeY; i++)
		for (j = 0; j < sizeX; j++)
			newPicture1[i + 2][j + 2] = (float)picture1[i * sizeX + j];
	float** newPicture2 = (float**)malloc(sizeof(float*) * (4 + sizeY));
	for (i = 0; i < sizeY + 4; i++)
		newPicture2[i] = (float*)malloc(sizeof(float) * (sizeX + 4));
	for (i = 0; i < sizeX + 4; i++)
		newPicture2[0][i] = newPicture2[1][i] = newPicture2[sizeY + 2][i] = newPicture2[sizeY + 3][i] = -1;
	for (i = 0; i < sizeY + 4; i++)
		newPicture2[i][0] = newPicture2[i][1] = newPicture2[i][sizeX + 2] = newPicture2[i][sizeX + 3] = -1;
	for (i = 0; i < sizeY; i++)
		for (j = 0; j < sizeX; j++)
			newPicture2[i + 2][j + 2] = (float)picture2[i * sizeX + j];
	for (i = 0; i < sizeY; i++)
	{
		for (j = 0; j < sizeX; j++)
		{
			float field1[25], field2[25];
			float avg1 = 0, avg2 = 0;
			int num = 0;
			for (k = i; k < i + 5; k++)
			{
				for (l = j; l < j + 5; l++)
				{
					if (newPicture1[k][l] != -1)
					{
						field1[num] = newPicture1[k][l];
						field2[num] = newPicture2[k][l];
						num++;
					}
				}
			}
			for (k = 0; k < num; k++)
			{
				avg1 += field1[k];
				avg2 += field2[k];
			}
			avg1 /= num;
			avg2 /= num;
			float temp1 = 0, temp2 = 0, temp3 = 0;
			for (k = 0; k < num; k++)
				temp1 += (field1[k] - avg1) * (field2[k] - avg2);
			for (k = 0; k < num; k++)
				temp2 += (field1[k] - avg1) * (field1[k] - avg1);
			for (k = 0; k < num; k++)
				temp3 += (field2[k] - avg2) * (field2[k] - avg2);
			Ncc[i][j] = temp1 / sqrt(temp2 * temp3);
			PD[i * sizeX + j] = (float)((1.0 - Ncc[i][j]) / 2.0);
			Ncc[i][j] = (float)((1.0 - Ncc[i][j]) / 2.0);
		}
	}
	for (i = 0; i < sizeY + 4; i++)
	{
		free(newPicture1[i]);
		free(newPicture2[i]);
	}
	free(newPicture1);
	free(newPicture2);
	qsort(PD, sizeX * sizeY, sizeof(PD[0]), cmp);
	float TH = PD[sizeY * sizeX / 2];
	for (k = 0; k < SuperpixelNum; k++)
	{
		int pictureX, pictureY;
		spoint* tempspspoint;
		tempspspoint = supPixel[k]->pixel;
		thisPixelNum = supPixel[k]->PixelNum;
		int bigNum = 0;
		for (i = 0; i < thisPixelNum; i++)
		{
			spoint tempPoint = tempspspoint[i];
			pictureX = tempPoint.x;
			pictureY = tempPoint.y;
			if (Ncc[pictureX][pictureY] > TH) bigNum++;
		}
		supPixel[k]->RegionChangeRate = (float)bigNum / thisPixelNum;
	}
	for (i = 0; i < sizeY; i++)
		free(Ncc[i]);
	free(Ncc);
	free(PD);
	float** Gradient = (float**)malloc(sizeof(float*) * sizeY);
	for (i = 0; i < sizeY; i++)
		Gradient[i] = (float*)malloc(sizeof(float) * sizeX);
	GetGradient(picture1, picture2, Gradient, sizeX, sizeY, SuperpixelNum, supPixel);
	for (k = 0; k < SuperpixelNum; k++)
	{
		float tempsum = 0;
		int pictureX, pictureY;
		spoint* tempspspoint;
		tempspspoint = supPixel[k]->pixel;
		for (i = 0; i < supPixel[k]->PixelNum; i++)
		{
			spoint tempPoint = tempspspoint[i];
			pictureX = tempPoint.x;
			pictureY = tempPoint.y;
			tempsum += Gradient[pictureX][pictureY] * Gradient[pictureX][pictureY];
		}
		supPixel[k]->SuperPixelDifferece = sqrt(tempsum);
	}
	for (i = 0; i < sizeY; i++)
		free(Gradient[i]);
	free(Gradient);
}

// function	release_superpixel
void release_superpixel(int SuperpixelNum, SuperPixel** supPixel)
{
	int i;
	for (i = 0; i < SuperpixelNum; i++)
	{
		free(supPixel[i]->pixel);
		free(supPixel[i]);
	}
	free(supPixel);
}

// function	BuildCostWeightMatrix
void BuildCostWeightMatrix(int** RegionAdjacencyMatrix, float** CostWeightMatrix, int MaxLabel, SuperPixel** supPixel, float weight1, float weight2, float weight3)
{
	for (int i = 0; i < MaxLabel; i++)
	{
		for (int j = 0; j < MaxLabel; j++)
		{
			if (RegionAdjacencyMatrix[i][j] == 1)
			{
				CostWeightMatrix[i][j] = supPixel[i]->RegionChangeRate + supPixel[j]->RegionChangeRate + weight1 * (supPixel[i]->Avg + supPixel[j]->Avg) + weight2 * fabs(supPixel[i]->Avg - supPixel[j]->Avg) + weight3 * (supPixel[i]->SuperPixelDifferece + supPixel[j]->SuperPixelDifferece);
			}
			else CostWeightMatrix[i][j] = 0;
		}
	}
}

// function	init
void init(queue<int>& q, int* Path)
{
	while (!q.empty())
		q.pop();
	for (int i = 0; i < MAXN; i++)
		Path[i] = -1;
}

// function	Min
float Min(float a, float b)
{
	return a < b ? a : b;
}

// function	bfs
float bfs(queue<int>& q, float** Net, int* Path, float* Lv, int& superpixel_num, int& Start, int& End)
{
	init(q, Path);
	Path[Start] = 0;
	Lv[Start] = (float)lmax;
	q.push(Start);
	while (!q.empty())
	{
		long t = q.front();
		q.pop();
		if (t == End) break;
		long i;
		for (i = 1; i <= superpixel_num; ++i)
		{
			if (i != Start && Path[i] == -1 && Net[t][i])
			{
				Lv[i] = Min(Lv[t], Net[t][i]);
				q.push(i);
				Path[i] = t;
			}
		}
	}
	if (Path[End] == -1) return -1;
	return Lv[End];
}

// function	cutpart
void cutpart(float** Net, int* Path, int& pre_superpixelNum, int* out, int& t)
{
	int i, j, k = 0;
	for (i = 0; i < MAXN; i++)
		for (j = 0; j < 2; j++)
			out[i * 2 + j] = -1;
	for (i = 1; i <= pre_superpixelNum; i++)
	{
		for (j = 1; j <= pre_superpixelNum; j++)
			if ((Net[j][i] || Net[i][j]) && Path[i] != -1 && Path[j] == -1)
			{
				out[k * 2] = i; out[k * 2 + 1] = j; k++;
			}
	}
}

// function	Ford_Fulkerson
void Ford_Fulkerson(queue<int>& q, float** Net, float** Net1, int* Path, float* Lv, int& superpixel_num, int& nn, int& Start, int& End, int& count1, int& count2, int* c1, int* c2, int& pre_superpixelNum, int* out, int& t, float** data)
{
	long i, j;
	nn = 0;
	for (i = 0; i < MAXN; i++)
		for (j = 0; j < MAXN; j++)
			Net[i][j] = Net1[i][j] = 0;
	float Max_Flow = 0;
	for (i = 1; i <= superpixel_num; ++i)
		for (j = 1; j <= superpixel_num; j++)
		{
			float cost;
			cost = data[i - 1][j - 1];
			if (cost != 0 && i <= j)
			{
				Net[i][j] = cost; nn++;
			}
			else Net1[i][j] = cost;
		}
	for (i = 0; i < count1; i++)
		Net[0][(int)c1[i]] = (float)99999999.0;
	for (i = 0; i < count2; i++)
		Net[(int)c2[i]][superpixel_num + 1] = (float)99999999.0;
	superpixel_num += 2;
	int kk = 0;
	for (i = 0; i < superpixel_num; i++)
		for (j = 0; j < superpixel_num; j++)
		{
			if (Net1[i][j] != 0 && j < i)
			{
				Net[i][superpixel_num + kk] = Net1[i][j];
				Net[superpixel_num + kk][j] = Net1[i][j];
				kk++;
			}
		}
	superpixel_num += nn;
	nn = nn * 3 + count1 + count2;
	Start = 0; End = pre_superpixelNum + 1;
	float step;
	while ((step = bfs(q, Net, Path, Lv, superpixel_num, Start, End)) != -1)//???????????
	{
		Max_Flow += step;
		int now = End;
		while (now != Start)
		{
			long pre = Path[now];
			Net[pre][now] -= step;
			Net[now][pre] += step;
			now = pre;
		}
	}
	cutpart(Net, Path, pre_superpixelNum, out, t);
}

// function	mincut
void mincut(int superpixel_num, int count1, int count2, float** data, int* c1, int* c2, int* out)
{
	int i;
	for (i = 0; i < count1; i++) c1[i] = c1[i] + 1;
	for (i = 0; i < count2; i++) c2[i] = c2[i] + 1;
	int nn;//edge
	int Start, End;
	int pre_superpixelNum = superpixel_num;
	int t = 100 * superpixel_num;
	queue<int> q;
	int* Path = (int*)malloc(sizeof(int) * MAXN);//Augmentation path
	float* Lv = (float*)malloc(sizeof(float) * t);//Augmented path capacity
	float** Net = (float**)malloc(sizeof(float*) * MAXN);//Residual network
	for (i = 0; i < MAXN; i++)
		Net[i] = (float*)malloc(sizeof(float) * MAXN);
	float** Net1 = (float**)malloc(sizeof(float*) * MAXN);//Residual network
	for (i = 0; i < MAXN; i++)
		Net1[i] = (float*)malloc(sizeof(float) * MAXN);
	Ford_Fulkerson(q, Net, Net1, Path, Lv, superpixel_num, nn, Start, End, count1, count2, c1, c2, pre_superpixelNum, out, t, data);
	free(Lv);
	free(Path);
	for (i = 0; i < MAXN; i++)
		free(Net[i]);
	free(Net);
	for (i = 0; i < MAXN; i++)
		free(Net1[i]);
	free(Net1);
}

// function	dfs_class
void dfs_class(int* cutline, int* Label, int* symbol, int* Labelclass, int curpos, int iWidth, int iHeight)
{
	int direction[4] = { -1,1,-iWidth,iWidth };
	stack<int> q;
	q.push(curpos);
	while (!q.empty())
	{
		int temp = q.top();
		q.pop();
		for (int k = 0; k < 4; k++)
		{
			if ((temp + direction[k]) < iWidth * iHeight && (temp + direction[k]) >= 0 && cutline[temp + direction[k]] == 0 && Label[temp + direction[k]] >= 0 && symbol[temp + direction[k]] == 0)
			{
				symbol[temp + direction[k]] = 1;
				Labelclass[temp + direction[k]] = 2;
				q.push(temp + direction[k]);
			}
		}
	}
}

std::pair<int*, int*> seamcut(int iWidth, int iHeight, int* mask, int* Label, int MaxLabel, std::vector<float> cost, int ImgsizeX1, int ImgsizeY1, int ImgsizeX2, int ImgsizeY2)
{
	//input the size of the image
	/*int iWidth = 472;
	int iHeight = 432;*/
	int WinSize = 20;
	double lambda = 0.1;
	//int* mask = new int[iWidth * iHeight];

	//ifstream f3("mask.txt");
	//for (int i = 0; i < iWidth * iHeight; i++)
	//{
	//	f3 >> mask[i];
	//}
	//f3.close();

	//int* Label = new int[iWidth * iHeight];

	//ifstream fl("label.txt");
	//for (int i = 0; i < iWidth * iHeight; i++)
	//{
	//	fl >> Label[i];
	//}
	//fl.close();

	//int MaxLabel = 262;

	float** CostWeightMatrix = new float* [MaxLabel];
	for (int k = 0; k < MaxLabel; k++)
		CostWeightMatrix[k] = new float[MaxLabel];

	//ifstream fc("cost.txt");
	//for (int i = 0; i < MaxLabel; i++)
	//	for (int j = 0; j < MaxLabel; j++)
	//	{
	//		fc >> CostWeightMatrix[i][j];
	//	}
	//fc.close();
	for (int i = 0; i < MaxLabel; i++) {
		for (int j = 0; j < MaxLabel; j++) {
			CostWeightMatrix[i][j] = cost[i * MaxLabel + j];
		}
	}

	int* classmask = new int[iWidth * iHeight];
	memset(classmask, 0, sizeof(int) * iWidth * iHeight);
	int num1 = 0, num2 = 0;

	//Intersection point
	//int ImgsizeX1 = 452, ImgsizeY1 = 84;
	//int ImgsizeX2 = 16, ImgsizeY2 = 278;

	classify(Label, classmask, ImgsizeX1, ImgsizeY1, ImgsizeX2, ImgsizeY2, iWidth, iHeight, num1, num2);
	int* class11 = new int[num1 + 1];
	memset(class11, 0, sizeof(int) * (num1 + 1));
	int* class22 = new int[num2 + 1];
	memset(class22, 0, sizeof(int) * (num2 + 1));
	int m = 0, n = 0;
	int max1 = 0, max2 = 0;
	for (int i = 0; i < iWidth * iHeight; i++)
	{
		if (classmask[i] == 1)
		{
			class11[m++] = Label[i];
			if (Label[i] > max1)
			{
				max1 = Label[i];
			}
		}
		else if (classmask[i] == 2)
		{
			class22[n++] = Label[i];
			if (Label[i] > max2)
			{
				max2 = Label[i];
			}
		}
		else continue;
	}
	int* temp1 = new int[max1 + 5];
	memset(temp1, 0, sizeof(int) * (max1 + 5));
	int* temp2 = new int[max2 + 5];
	memset(temp2, 0, sizeof(int) * (max2 + 5));
	for (int i = 0; i < m; i++)
	{
		temp1[class11[i]]++;
	}
	for (int i = 0; i < n; i++)
	{
		temp2[class22[i]]++;
	}
	if (Label[ImgsizeY1 * iWidth + ImgsizeX1] >= 0 && Label[ImgsizeY1 * iWidth + ImgsizeX1] < max1 + 1)
	{
		temp1[Label[ImgsizeY1 * iWidth + ImgsizeX1]] = 0;
		temp2[Label[ImgsizeY1 * iWidth + ImgsizeX1]] = 0;
	}
	if (Label[ImgsizeY2 * iWidth + ImgsizeX2] >= 0 && Label[ImgsizeY2 * iWidth + ImgsizeX2] < max2 + 1)
	{
		temp1[Label[ImgsizeY2 * iWidth + ImgsizeX2]] = 0;
		temp2[Label[ImgsizeY2 * iWidth + ImgsizeX2]] = 0;
	}
	int class1num = 0, class2num = 0;
	for (int i = 0; i < max1 + 1; i++)
	{
		if (temp1[i] != 0)
		{
			class1num++;
		}
	}
	for (int i = 0; i < max2 + 1; i++)
	{
		if (temp2[i] != 0)
		{
			class2num++;
		}
	}
	int* class1 = new int[class1num + 1];
	int* class2 = new int[class2num + 1];
	int x = 0, y = 0;
	for (int i = 0; i < max1 + 1; i++)
	{
		if (temp1[i] != 0)
		{
			class1[x++] = i;
		}
	}
	for (int i = 0; i < max2 + 1; i++)
	{
		if (temp2[i] != 0)
		{
			class2[y++] = i;
		}
	}

	int* out = (int*)malloc(sizeof(int) * MAXN * 2);
	//cutline
	int* cutline = new int[iWidth * iHeight];
	int dd[4] = { -1,-iWidth,1,iWidth };
	memset(cutline, 0, sizeof(int) * iWidth * iHeight);

	double dur;
	clock_t start, end;
	start = clock();
	mincut(MaxLabel, class1num, class2num, CostWeightMatrix, class1, class2, out);
	end = clock();
	dur = (double)(end - start);
	//cout << dur << endl;
	//cout << dur / CLOCKS_PER_SEC << endl;

	for (int i = 0; i < 2 * MaxLabel && out[i * 2] != -1; i++)
	{
		for (int j = 0; j < iWidth * iHeight; j++)
		{
			for (int k = 0; k < 4; k++)
			{	
				int abs_idx_labels = abs(j + dd[k]);
				if (Label[j] == (out[i * 2] - 1) && Label[abs_idx_labels] == (out[i * 2 + 1] - 1)) {
					cutline[j] = 1;
				}
			}
		}
	}

	int dd2[8] = { -iWidth - 1,-1,iWidth - 1,iWidth,iWidth + 1,1,-iWidth + 1,-iWidth };
	int* tag1 = new int[iWidth * iHeight];
	memset(tag1, 0, sizeof(int) * iWidth * iHeight);
	int* tag2 = new int[iWidth * iHeight];
	memset(tag2, 0, sizeof(int) * iWidth * iHeight);
	int spos1, spos2;
	int scount1 = 0, scount2 = 0;
	for (int i = 0; i < iWidth * iHeight; i++)
	{
		if (cutline[i] == 1)
		{
			spos1 = i;
			tag1[spos1] = 1;
			break;
		}
	}
	for (int i = 0; i < iWidth * iHeight; i++)
	{
		if (cutline[iWidth * iHeight - 1 - i] == 1)
		{
			spos2 = iWidth * iHeight - 1 - i;
			tag2[spos2] = 1;
			break;
		}
	}
	stack<int> q1, q2;
	q1.push(spos1);
	while (!q1.empty())
	{
		int stemp1 = q1.top();
		q1.pop();
		for (int k = 0; k < 8; k++)
		{
			if ((stemp1 + dd2[k]) < iWidth * iHeight && (stemp1 + dd2[k]) >= 0 && cutline[stemp1 + dd2[k]] == 1 && tag1[stemp1 + dd2[k]] == 0)
			{
				tag1[stemp1 + dd2[k]] = 1;
				scount1++;
				q1.push(stemp1 + dd2[k]);
			}
		}
	}
	q2.push(spos2);
	while (!q2.empty())
	{
		int stemp2 = q2.top();
		q2.pop();
		for (int k = 0; k < 8; k++)
		{
			if ((stemp2 + dd2[k]) < iWidth * iHeight && (stemp2 + dd2[k]) >= 0 && cutline[stemp2 + dd2[k]] == 1 && tag2[stemp2 + dd2[k]] == 0)
			{
				tag2[stemp2 + dd2[k]] = 1;
				scount2++;
				q2.push(stemp2 + dd2[k]);
			}
		}
	}
	for (int i = 0; i < iWidth * iHeight; i++)
	{
		if (scount1 > scount2)
		{
			if (tag1[i] == 0)
				cutline[i] = 0;
		}
		else
		{
			if (tag2[i] == 0)
				cutline[i] = 0;
		}
	}

	//ofstream xian("cutline.txt");
	//for (int i = 0; i < iHeight; i++)
	//{
	//	for (int j = 0; j < iWidth; j++)
	//	{
	//		xian << cutline[i * iWidth + j] << " ";
	//	}
	//	xian << endl;
	//}
	//xian.close();


	int* Labelclass = new int[iWidth * iHeight];
	memset(Labelclass, 0, sizeof(int) * iWidth * iHeight);
	for (int i = 0; i < iWidth * iHeight; i++)
	{
		Labelclass[i] = 0;
	}
	int* symbol = new int[iWidth * iHeight];
	memset(symbol, 0, sizeof(int) * iWidth * iHeight);
	int curpos;
	for (int i = 0; i < iWidth * iHeight; i++)
	{
		if (cutline[i] == 0 && Label[i] >= 0 && symbol[i] == 0)
		{
			curpos = i;
			symbol[curpos] = 1;
			Labelclass[curpos] = 2;
			break;
		}
	}
	dfs_class(cutline, Label, symbol, Labelclass, curpos, iWidth, iHeight);
	for (int i = 0; i < iWidth * iHeight; i++)
	{
		if (Label[i] >= 0 && symbol[i] == 0)
		{
			Labelclass[i] = 1;
			symbol[i] = 1;
		}
	}

	for (int i = 1; i < iWidth - 1; i++)
	{
		if (mask[i] == 1)
			Labelclass[i] = Labelclass[iWidth + i];
		if (mask[(iHeight - 1) * iWidth + i] == 1)
			Labelclass[(iHeight - 1) * iWidth + i] = Labelclass[(iHeight - 2) * iWidth + i];
	}
	for (int j = 1; j < iHeight - 1; j++)
	{
		if (mask[j * iWidth] == 1)
			Labelclass[j * iWidth] = Labelclass[j * iWidth + 1];
		if (mask[j * iWidth + iWidth - 1] == 1)
			Labelclass[j * iWidth + iWidth - 1] = Labelclass[j * iWidth + iWidth - 2];
	}
	if (mask[0] == 1)
		Labelclass[0] = Labelclass[iWidth + 1];
	if (mask[iWidth - 1] == 1)
		Labelclass[iWidth - 1] = Labelclass[2 * (iWidth - 1)];
	if (mask[(iHeight - 1) * iWidth] == 1)
		Labelclass[(iHeight - 1) * iWidth] = Labelclass[(iHeight - 2) * iWidth + 1];
	if (mask[iHeight * iWidth - 1] == 1)
		Labelclass[iHeight * iWidth - 1] = Labelclass[(iHeight - 1) * iWidth - 2];

	//ofstream outlabel("labelclass.txt");
	//for (int i = 0; i < iHeight; i++)
	//{
	//	for (int j = 0; j < iWidth; j++)
	//	{
	//		outlabel << Labelclass[i * iWidth + j] << " ";
	//	}
	//	outlabel << endl;
	//}
	//outlabel.close();

	delete[]classmask;
	delete[]class11;
	delete[]class22;
	delete[]temp1;
	delete[]temp2;
	delete[]class1;
	delete[]class2;
	for (int i = 0; i < MaxLabel; i++)
	{
		delete[]CostWeightMatrix[i];
	}
	delete[]CostWeightMatrix;
	delete[]Label; Label = NULL;
	free(out);
	//delete[]cutline;
	delete[]symbol;
	delete[]tag1;
	delete[]tag2;

	pair<int*, int*> result;
	result = make_pair(cutline, Labelclass);
	return result;
}