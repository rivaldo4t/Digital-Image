#pragma once
#include <iostream>
#include <fstream>
#include <GL/glut.h>
#include <memory>
#include <vector>
#include <numeric>
#include <random>
#include <chrono>
#include "utility.hpp"
#include "image.hpp"
#include "ppm.hpp"
#define _USE_MATH_DEFINES
#include <math.h>
#include "matrix.hpp"

// #define CONVOLUTION_FILTERS
// #define RASTERIZED_SHAPES
// #define PROCEDURAL_IMAGE_GEN
// #define TRANSFORMATIONS
// #define COMPOSITIONS

#define DITHERING
// #define ANTIALIASING

int width, height;
std::vector<uint8_t> pixmap;
std::vector<uint8_t> pixmap2;
std::vector<uint8_t> pixmap3;

void clampImgBound(int& row, int& col)
{
	// mirror
	if (row < 0)
		row = 0;
	else if (row > height - 1)
		row = height - 1;
	
	if (col < 0)
		col = 0;
	else if (col > width - 1)
		col = width - 1;
}

class Kernel
{
public:
	int kHeight, kWidth;
	std::vector<float> weights;
	Kernel(int sizeX, int sizeY) : kHeight(sizeX), kWidth(sizeY) {} // keep sizeX, sizeY odd for now
};

class ConvolutionFilter : public Kernel
{
public:
	// hold partition of unity
	ConvolutionFilter(int sizeX, int sizeY) : Kernel(sizeX, sizeY) {}

	void boxFilter()
	{
		weights.resize(kHeight * kWidth, 1.0 / float(kHeight * kWidth));
	}

	void radialFilter()
	{
		weights.resize(kHeight * kWidth, 0.0);

		float sum = 0;
		for (int i = 0; i < kHeight; ++i)
		{
			for (int j = 0; j < kWidth; ++j)
			{
				int index = i * kWidth + j;
				int centerX = kHeight / 2;
				int centerY = kWidth / 2;
				float radius = sqrt(float(centerX * centerX + centerY * centerY));
				float dist = sqrt(float(pow(centerX - i, 2.0) + pow(centerY - j, 2.0)));

				weights[index] = pow(radius - dist, 2);
				sum += weights[index];
			}
		}

		for (auto& k : weights)
			k /= sum;
	}

	void motionFilter()
	{
		weights.resize(kHeight * kWidth, 0.0);
		float sum = 0;

		for (int i = 0; i < kHeight; ++i)
		{
			for (int j = 0; j < kWidth; ++j)
			{
				int index = i * kWidth + j;

				// enable any one of the below

				// horizontal
				// weights[index] = pow(j - kWidth / 2, 2) + kWidth / 2;

				// vertical
				// weights[index] = pow(j - kHeight / 2, 2) + kHeight / 2;

				// diagonal
				if (i == j)
					weights[index] = kHeight;

				sum += abs(weights[index]);
			}
		}

		for (auto& k : weights)
			k /= sum;
	}

	Color eval(int x, int y, int op = 0)
	{
		Color c;
		
		for (int i = 0; i < kHeight; ++i)
		{
			for (int j = 0; j < kWidth; ++j)
			{
				int row = x + i - 0.5 * (kHeight - 1);
				int col = y + j - 0.5 * (kWidth - 1);
				clampImgBound(row, col);
				int pix = (row * width + col) * 3;
				Color neighborColor(pixmap[pix + 0] / 255.0, pixmap[pix + 1] / 255.0, pixmap[pix + 2] / 255.0);
				
				if (op == 1)
					neighborColor = Color(pixmap2[pix + 0] / 255.0, pixmap2[pix + 1] / 255.0, pixmap2[pix + 2] / 255.0);

				int index = i * kWidth + j;

				c += neighborColor * weights[index];
			}
		}

		return c;
	}
};

class MorphologicalFilter : public Kernel
{
	bool isDilation;
	float maxWeight, minWeight;
public:
	MorphologicalFilter(int sizeX, int sizeY) : Kernel(sizeX, sizeY) {}

	void dilationFilter()
	{
		isDilation = true;
		weights.resize(kHeight * kWidth, 0.0);
		maxWeight = INT_MIN;

		for (int i = 0; i < kHeight; ++i)
		{
			for (int j = 0; j < kWidth; ++j)
			{
				int index = i * kWidth + j;
				int centerX = kHeight / 2;
				int centerY = kWidth / 2;
				float radius = sqrt(float(centerX * centerX + centerY * centerY));
				float dist = sqrt(float(pow(centerX - i, 2.0) + pow(centerY - j, 2.0)));

				//weights[index] = abs(radius - dist) + kHeight / 2;
				weights[index] = 1;
				maxWeight = std::max(maxWeight, weights[index]);
			}
		}
	}

	void erosionFilter()
	{
		weights.resize(kHeight * kWidth, 0.0);
		isDilation = false;
		minWeight = INT_MAX;

		for (int i = 0; i < kHeight; ++i)
		{
			for (int j = 0; j < kWidth; ++j)
			{
				int index = i * kWidth + j;
				int centerX = kHeight / 2;
				int centerY = kWidth / 2;
				float radius = sqrt(float(centerX * centerX + centerY * centerY));
				float dist = sqrt(float(pow(centerX - i, 2.0) + pow(centerY - j, 2.0)));

				//weights[index] = abs(radius - dist) + kHeight / 2;
				weights[index] = 1;
				minWeight = std::min(minWeight, weights[index]);
			}
		}
	}

	Color eval(int x, int y)
	{
		Color c;
		if (!isDilation)
			c = Color(1.0, 1.0, 1.0);

		for (int i = 0; i < kHeight; ++i)
		{
			for (int j = 0; j < kWidth; ++j)
			{
				int row = x + i - 0.5 * (kHeight - 1);
				int col = y + j - 0.5 * (kWidth - 1);
				clampImgBound(row, col);
				int pix = (row * width + col) * 3;
				Color neighborColor(pixmap[pix + 0] / 255.0, pixmap[pix + 1] / 255.0, pixmap[pix + 2] / 255.0);

				int index = i * kWidth + j;
				if (neighborColor.r > 0)
				{
					index = index;
				}
				Color c2 = neighborColor * weights[index];
				if (isDilation)
					c = c2 > c ? c2 : c;
				else
					c = c2 < c ? c2 : c;
			}
		}

		if (isDilation)
			c = c / maxWeight;
		else 
			c = c / minWeight;
		return c;
	}
};

class DerivativeFilter : public Kernel
{
	float absWeightSum;
	std::vector<float> weights2;
	float absWeightSum2;
	bool isEdge = false;
public:
	DerivativeFilter(int sizeX, int sizeY) : Kernel(sizeX, sizeY) {}
	
	void embossFilter1()
	{
		weights.resize(kHeight * kWidth, 0.0);
		absWeightSum = 0;
		float sum = 0;

		for (int i = 0; i < kHeight; ++i)
		{
			for (int j = 0; j < kWidth; ++j)
			{
				int index = i * kWidth + j;

				if (j <= kWidth / 2)
					weights[index] = j + 1;

				sum += weights[index];
			}
		}

		for (auto& k : weights)
		{
			k -= sum / float(kWidth * kHeight);
			absWeightSum += abs(k);
		}
	}

	void embossFilter2()
	{
		weights.resize(kHeight * kWidth, 0.0);

		weights[0] = 1;
		weights[kWidth - 1] = -1;
		weights[kWidth * (kHeight - 1)] = -1;
		weights[kWidth * (kHeight - 1) + kWidth - 1] = 1;

		absWeightSum = 4;
	}

	void embossFilter3()
	{
		weights.resize(kHeight * kWidth, 0.0);
		absWeightSum = 0;

		for (int i = 0; i < kHeight; ++i)
		{
			for (int j = 0; j < kWidth; ++j)
			{
				int index = i * kWidth + j;

				if (i == kHeight / 2 && j == kWidth / 2)
					weights[index] = -(kHeight * kWidth - 1);
				else
					weights[index] = 1;

				absWeightSum += abs(weights[index]);
			}
		}
	}

	void embossFilter4()
	{
		/*weights = { 0, -1, 0, -1, 5, -1, 0, -1, 0 };
		absWeightSum = 9;
		return;*/
		weights.resize(kHeight * kWidth, 0.0);
		int i, j;
		i = 0, j = 1;
		weights[i * kWidth + j] = 1;
		i = 1, j = 0;
		weights[i * kWidth + j] = 1;
		i = 1, j = 1;
		weights[i * kWidth + j] = -4;
		i = 2, j = 2;
		weights[i * kWidth + j] = 1;
		i = 1, j = 2;
		weights[i * kWidth + j] = 1;
		absWeightSum = 8;
	}

	// try thresholding for better results
	void edgeFilter()
	{
		isEdge = true;
		
		weights = { -1, -2, -1, 0, 0, 0, 1, 2, 1 };
		absWeightSum = 8;
		weights2 = { 1, 0, -1, 2, 0, -2, 1, 0, -1 };
		absWeightSum2 = 8;
		return;

		// for Gradient in Y
		weights.resize(kHeight * kWidth, 0.0);
		absWeightSum = 0;

		for (int i = 0; i < kHeight; ++i)
		{
			for (int j = 0; j < kWidth; ++j)
			{
				int index = i * kWidth + j;

				weights[index] = kHeight / 2 - i;

				absWeightSum += abs(weights[index]);
			}
		}

		// for Gradient in X
		weights2.resize(kHeight * kWidth, 0.0);
		absWeightSum2 = 0;

		for (int j = 0; j < kWidth; ++j)
		{
			for (int i = 0; i < kHeight; ++i)
			{
				int index = i * kWidth + j;

				weights2[index] = kWidth/2 - j;

				absWeightSum2 += abs(weights2[index]);
			}
		}
	}

	Color evalEdgeFilter(int x, int y)
	{
		Color cGradY, cGradX, cGradMag;

		for (int i = 0; i < kHeight; ++i)
		{
			for (int j = 0; j < kWidth; ++j)
			{
				int row = x + i - 0.5 * (kHeight - 1);
				int col = y + j - 0.5 * (kWidth - 1);
				clampImgBound(row, col);
				int pix = (row * width + col) * 3;
				Color neighborColor(pixmap[pix + 0] / 255.0, pixmap[pix + 1] / 255.0, pixmap[pix + 2] / 255.0);

				int index = i * kWidth + j;

				cGradY += neighborColor * weights[index];
				cGradX += neighborColor * weights2[index];
			}
		}
		
#if 0
		double tempR = (abs(cGradY.r) + abs(cGradX.r)) / 0.5;
		tempR = 1 - tempR;
		double tempG = (abs(cGradY.g) + abs(cGradX.g)) / 0.5;
		tempG = 1 - tempG;
		double tempB = (abs(cGradY.b) + abs(cGradX.b)) / 0.5;
		tempB = 1 - tempB;
		return Color(tempR, tempR, tempR);
#endif

		// you do this to get only the direction; I guess!
#if 0
		cGradY = Color(	(cGradY.r + 0.5*absWeightSum) / absWeightSum,
						(cGradY.g + 0.5*absWeightSum) / absWeightSum,
						(cGradY.b + 0.5*absWeightSum) / absWeightSum);

		cGradX = Color(	(cGradX.r + 0.5*absWeightSum2) / absWeightSum2,
						(cGradX.g + 0.5*absWeightSum2) / absWeightSum2,
						(cGradX.b + 0.5*absWeightSum2) / absWeightSum2);
#endif

#if 0
		cGradMag = cGradY;
		cGradMag += cGradX;
		cGradMag = cGradMag * 0.5;
#else

		cGradMag = Color(sqrt(cGradY.r * cGradY.r + cGradX.r * cGradX.r), 
			sqrt(cGradY.g * cGradY.g + cGradX.g * cGradX.g),
			sqrt(cGradY.b * cGradY.b + cGradX.b * cGradX.b));
#endif
		cGradMag = Color(1 - cGradMag.r, 1 - cGradMag.r, 1 - cGradMag.r);


		return cGradMag;
	}

	Color eval(int x, int y)
	{
		if (isEdge)
			return evalEdgeFilter(x, y);

		Color c;
		for (int i = 0; i < kHeight; ++i)
		{
			for (int j = 0; j < kWidth; ++j)
			{
				int row = x + i - 0.5 * (kHeight - 1);
				int col = y + j - 0.5 * (kWidth - 1);
				clampImgBound(row, col);
				int pix = (row * width + col) * 3;
				Color neighborColor(pixmap[pix + 0] / 255.0, pixmap[pix + 1] / 255.0, pixmap[pix + 2] / 255.0);

				int index = i * kWidth + j;

				c += neighborColor * weights[index];
			}
		}
		
		return Color(	(c.r + 0.5*absWeightSum) / absWeightSum,
						(c.g + 0.5*absWeightSum) / absWeightSum,
						(c.b + 0.5*absWeightSum) / absWeightSum		);
	}
};

double pieceWiseLinearInterpolation(std::vector<double>& curveParam, double c)
{
// stored as (r0, r1, r2, r3 ..., rn)
// points are given as (0, r0), (1/n, r1), (2/n, r2), (3/n, r3) and (1, rn)
	size_t n = curveParam.size();
	if (n <= 1)
		return c;

	double delta_x = 1 / double(n-1);
	int i;
	for (i = 0; i < n; ++i)
	{
		if (c < i * delta_x)
			break;
	}
	if (i >= n - 1)
		i = n - 2;

	double xk = i * delta_x;
	double xk_1 = (i + 1)*delta_x;
	double yk = curveParam[i];
	double yk_1 = curveParam[i + 1];

	double interpolated_c = ((xk_1 - c) * yk + (c - xk) * yk_1) / delta_x;
	return interpolated_c;
}

void curvesAdjustment(Image& img)
{
	width = img.w;
	height = img.h;
	img.flip();
	pixmap = std::move(img.data);

	// ----------------------------------------- curves adjustments -----------------------------------------

	// reduce midtones
	std::vector<double> curveParam   = { 0.0, 0.4, 1.0 };
	
	// ig - crema filter
	std::vector<double> curveParam_r = { 0.0, 0.5, 1.0 };
	std::vector<double> curveParam_g = { 0.1, 0.5, 1.0 };
	std::vector<double> curveParam_b = { 0.0, 0.6, 1.0 };
	
	// increase contrast
	// note how 0.2 < 0.25 and 0.8 > 0.75 to make dark regions darker and bright regions brighter
	curveParam = {0.0, 0.2, 0.5, 0.8, 1.0};

	// another filter
	curveParam_r = { 0.1, 0.9 };
	curveParam_g = { 0.1, 0.6, 1.0 };
	curveParam_b = { 0.0, 0.6, 0.65, 0.7, 0.9 };

	// ------------------------------------------------------------------------------------------------------

#pragma omp parellel for
	for (int i = 0; i < height; ++i)
	{
		for (int j = 0; j < width; ++j)
		{
			int p = (i * width + j) * 3;

			double r = pixmap[p + 0] / 255.0;
			double g = pixmap[p + 1] / 255.0;
			double b = pixmap[p + 2] / 255.0;

			// either adjust the curve per channel or for all of them
			// use corresponding curveParam vectors to manipulate
			r = clamp<double>(pieceWiseLinearInterpolation(curveParam_r, r));
			g = clamp<double>(pieceWiseLinearInterpolation(curveParam_g, g));
			b = clamp<double>(pieceWiseLinearInterpolation(curveParam_b, b));
			
			pixmap[p + 0] = (uint8_t)(r * 255);
			pixmap[p + 1] = (uint8_t)(g * 255);
			pixmap[p + 2] = (uint8_t)(b * 255);
		}
	}
}

void replaceHues(Image& src, Image& dst)
{
	if (src.w != dst.w || src.h != dst.h)
	{
		std::cout << "Please make sure to match source and destination image dimensions\n";
		return;
	}

	width = src.w;
	height = src.h;
	src.flip();
	dst.flip();
	pixmap = dst.data;

#pragma omp parellel for
	for (int i = 0; i < height; ++i)
	{
		for (int j = 0; j < width; ++j)
		{
			int p = (i * width + j) * 3;

			double h_src, s_src, v_src;
			RGBtoHSV(src.data[p + 0], src.data[p + 1], src.data[p + 2], h_src, s_src, v_src);

			double h_dst, s_dst, v_dst;
			RGBtoHSV(dst.data[p + 0], dst.data[p + 1], dst.data[p + 2], h_dst, s_dst, v_dst);

			int r, g, b;
			HSVtoRGB(r, g, b, h_src, s_dst, v_dst);

			pixmap[p + 0] = (uint8_t)(r);
			pixmap[p + 1] = (uint8_t)(g);
			pixmap[p + 2] = (uint8_t)(b);
		}
	}
}

void inverseTransform(double& X, double& Y)
{
	Vector v(X, Y, 1);
	std::shared_ptr<Matrix> invM;
	bool wrap = false;
	
	// Rotation
	float theta = -30 * 3.14 / 180.0;
	//invM = std::make_shared<RotationMatrix>(theta);

	// Scaling
	float sx = 0.5, sy = 0.5;
	//invM = std::make_shared<ScaleMatrix>(sx, sy);
	
	// Mirror
	float mx = -1, my = 1;
	//invM = std::make_shared<ScaleMatrix>(mx, my);
	//wrap = true;

	// Translation
	float tx = -100, ty = -200;
	//invM = std::make_shared<TranslateMatrix>(tx, ty);

	// Shear
	float shx = 1, shy = 0;
	//invM = std::make_shared<ShearMatrix>(shx, shy);

	// Perspective
	float px = 0.000, py = 0.003;
	invM = std::make_shared<PerspectiveMatrix>(px, py);
	
	// Invert the matrix and apply transformation
	invM->invertMatrix();
	v = (*invM.get() * v);
	v.normalize();
	
	// scale down to see other transformations clearly
#if 0
	invM = std::make_shared<ScaleMatrix>(0.5, 0.5);
	invM->invertMatrix();
	v = (*invM.get() * v);
	v.normalize();
#endif

	// translate to see other transformations clearly
#if 0	
	invM = std::make_shared<TranslateMatrix>(400, 100);
	invM->invertMatrix();
	v = (*invM.get() * v);
	v.normalize();
#endif

	X = v[0];
	Y = v[1];
	// wrap edges for mirror
	if (wrap)
	{
		X = int(v[0] + width) % width;
		Y = int(v[1] + height) % height;
	}
}

void getPerspectiveTransformedCorners(std::vector<std::vector<double>>& transFormedcorners)
{
	Vector vec0(0, 0, 1), vec1(0, height - 1, 1), vec2(width - 1, height - 1, 1), vec3(width - 1, 0, 1);
	std::vector<Vector> cornerVectors{ vec0, vec1, vec2, vec3 };
	
	std::shared_ptr<Matrix> invM;

#if 1
	invM = std::make_shared<PerspectiveMatrix>(0.0004, 0.0004);
	for (auto& corner : cornerVectors)
	{
		corner = (*invM.get() * corner);
		corner.normalize();
	}
#endif
#if 0
	invM = std::make_shared<ScaleMatrix>(1, 1);
	for (auto& corner : cornerVectors)
	{
		corner = (*invM.get() * corner);
		corner.normalize();
	}
#endif
#if 0
	invM = std::make_shared<TranslateMatrix>(0, 0);
	for (auto& corner : cornerVectors)
	{
		corner = (*invM.get() * corner);
		corner.normalize();
	}
#endif

	transFormedcorners.push_back({cornerVectors[0][0], cornerVectors[1][0], cornerVectors[2][0], cornerVectors[3][0] });
	transFormedcorners.push_back({ cornerVectors[0][1], cornerVectors[1][1], cornerVectors[2][1], cornerVectors[3][1] });
}

void bilinearWarpTransform(double& X, double& Y, std::vector<std::vector<double>>& transFormedcorners)
{
	std::vector<double> x = transFormedcorners[0], y = transFormedcorners[1];

	double a0, a1, a2, a3, b0, b1, b2, b3;
	a0 = x[0];
	a1 = x[3] - x[0];
	a2 = x[1] - x[0];
	a3 = x[2] - x[1] - x[3] + x[0];
	b0 = y[0];
	b1 = y[3] - y[0];
	b2 = y[1] - y[0];
	b3 = y[2] - y[1] - y[3] + y[0];
	
	double c0, c1, c2;
	c0 = a1 * (b0 - Y) + b1 * (X - a0);
	c1 = a3 * (b0 - Y) + b3 * (X - a0) + a1 * b2 - a2 * b1;
	c2 = a3 * b2 - a2 * b3;
	
	double u1, v1, u2, v2, u, v;
	v1 = -0.5 * (c1 / c2 - sqrt(c1 * c1 - 4 * c2 * c0) / c2);
	v2 = -0.5 * (c1 / c2 + sqrt(c1 * c1 - 4 * c2 * c0) / c2);
	u1 = (X - a0 - a2 * v1) / (a1 + a3 * v1);
	u2 = (X - a0 - a2 * v2) / (a1 + a3 * v2);
	
	if (u1 >= 0 && u1 <= 1 && v1 >= 0 && v1 <= 1)
	{
		u = u1;
		v = v1;
	}
	else if (u2 >= 0 && u2 <= 1 && v2 >= 0 && v2 <= 1)
	{
		u = u2;
		v = v2;
	}
	else
	{
		X = -1;
		Y = -1;
		return;
	}
	
	X = (int)((width - 1)*u + 0.5);
	Y = (int)((height - 1)*v + 0.5);
}

void warpTransform(double& X, double& Y)
{
#if 0
	// simple sine inverse transform
	double x = X, y = Y;
	double amp = 16;
	X = x - amp * sin(y / amp);
	Y = y - amp * sin(x / amp);
#else
	// twirl inverse transform
	double a = 10, b = 50; // a = rotation of twirl, b = size of twirl
	double twirlX = 400, twirlY = 200; // location of twirl
	double x = X - twirlX, y = Y - twirlY;
	double angle = abs(a*exp(-(x*x + y*y) / (b*b)));
	double u = cos(angle)*x + sin(angle)*y;
	double v = -sin(angle)*x + cos(angle)*y;
	X = u + twirlX;
	Y = v + twirlY;
#endif
}

// init in main
ConvolutionFilter k(3, 3);
Color compositeOperation(int pix, int i = 0, int j = 0)
{
	// bg
	Color c0(pixmap[pix + 0] / 255.0, pixmap[pix + 1] / 255.0, pixmap[pix + 2] / 255.0);
	// fg
	Color c1(pixmap2[pix + 0] / 255.0, pixmap2[pix + 1] / 255.0, pixmap2[pix + 2] / 255.0);
	// alpha mask
	Color c2(pixmap3[pix + 0] / 255.0, pixmap3[pix + 1] / 255.0, pixmap3[pix + 2] / 255.0);

	double alpha0 = 0.5, alpha1 = 0.5;
	Color compose;
	int op = 8;
	
	// Add support for associative over operation for multiple layer blending
	
	switch (op)
	{
		// Normal		
		case 1: compose = c1; break;
		
		// Multiply
		case 2: compose = c0 * c1; break;
		
		// Darken
		case 3:compose = min(c0, c1); break;
		
		// Linear Dodge
		case 4:compose = c1 + c0; break;
		
		// Lighter Color
		case 5:compose = max(c0, c1); break;
		
		// Difference
		case 6:compose = c1 - c0; break;
		
		// Exclusion
		case 7:compose = Color(1.0) - c0 * c1; break;
		
		// Masking using Green Screen
		case 8: {
			// direct masking with bg and fg pixel color values
			// gives sharper masks
			//compose = c1; alpha1 = c2.r < 0.6 ? 0.0 : 1.0; break;
			
			// avg of bg and fg over a neighborhood of 3x3 or 9x9 pixels
			// using convolution filter for calculating avg
			// gives blurred mask with green halo
			compose = c1; alpha1 = c2.r;
			if (c2.r < 1.0 && c2.r > 0.0)
			{
				Color c00 = k.eval(i, j);
				Color c11 = k.eval(i, j, 1);
				compose = c11;
				c0 = c00;
				break;
				Color c = c11 * alpha0 + c00 * (1 - alpha0);
				Color a = (c - c00) / (c11 - c00);
				if (a.r < 1 && a.r > 0)
					alpha1 = a.r;
				else if (a.g < 1 && a.g > 0)
					alpha1 = a.g;
				else
					alpha1 = a.b;
			}
			break;
		}
		
		// No composition
		default:compose = c0; break;
	}

	compose.clamp();
	return compose * alpha1 + c0 * (1 - alpha1);
}

Color floydSteinbergDither(int i, int j)
{
	// diffuses the error to 4 unprocessed neighboring pixels

	int pix = (i * width + j) * 3;
	int ny, nx, nc, np;
	float oldColorR, oldColorG, oldColorB;
	float newColorR, newColorG, newColorB;
	float quant_errorR, quant_errorG, quant_errorB;

	oldColorR = pixmap[pix] / 255.0;
	if (oldColorR < 0.25)
		newColorR = 0;
	else if (oldColorR >= 0.25 && oldColorR < 0.5)
		newColorR = 0.3;
	else if (oldColorR >= 0.5 && oldColorR < 0.75)
		newColorR = 0.7;
	else
		newColorR = 1.0;

	oldColorG = pixmap[pix + 1] / 255.0;
	if (oldColorG < 0.25)
		newColorG = 0;
	else if (oldColorG >= 0.25 && oldColorG < 0.5)
		newColorG = 0.3;
	else if (oldColorG >= 0.5 && oldColorG < 0.75)
		newColorG = 0.7;
	else
		newColorG = 1.0;

	oldColorB = pixmap[pix + 2] / 255.0;
	if (oldColorB < 0.25)
		newColorB = 0;
	else if (oldColorB >= 0.25 && oldColorB < 0.5)
		newColorB = 0.3;
	else if (oldColorB >= 0.5 && oldColorB < 0.75)
		newColorB = 0.7;
	else
		newColorB = 1.0;

	quant_errorR = oldColorR - newColorR;
	quant_errorG = oldColorG - newColorG;
	quant_errorB = oldColorB - newColorB;

	ny = i; nx = j + 1;
	np = (ny * width + nx) * 3;
	nc = pixmap[np];
	nc += quant_errorR * 7 / 16 * 255;
	pixmap[np] = nc > 255 ? 255 : nc;
	nc = pixmap[np + 1];
	nc += quant_errorG * 7 / 16 * 255;
	pixmap[np + 1] = nc > 255 ? 255 : nc;
	nc = pixmap[np + 2];
	nc += quant_errorB * 7 / 16 * 255;
	pixmap[np + 2] = nc > 255 ? 255 : nc;

	ny = i + 1; nx = j - 1;
	np = (ny * width + nx) * 3;
	nc = pixmap[np];
	nc += quant_errorR * 3 / 16 * 255;
	pixmap[np] = nc > 255 ? 255 : nc;
	nc = pixmap[np + 1];
	nc += quant_errorG * 3 / 16 * 255;
	pixmap[np + 1] = nc > 255 ? 255 : nc;
	nc = pixmap[np + 2];
	nc += quant_errorB * 3 / 16 * 255;
	pixmap[np + 2] = nc > 255 ? 255 : nc;

	ny = i + 1; nx = j;
	np = (ny * width + nx) * 3;
	nc = pixmap[np];
	nc += quant_errorR * 5 / 16 * 255;
	pixmap[np] = nc > 255 ? 255 : nc;
	nc = pixmap[np + 1];
	nc += quant_errorG * 5 / 16 * 255;
	pixmap[np + 1] = nc > 255 ? 255 : nc;
	nc = pixmap[np + 2];
	nc += quant_errorB * 5 / 16 * 255;
	pixmap[np + 2] = nc > 255 ? 255 : nc;

	ny = i + 1; nx = j + 1;
	np = (ny * width + nx) * 3;
	nc = pixmap[np];
	nc += quant_errorR * 1 / 16 * 255;
	pixmap[np] = nc > 255 ? 255 : nc;
	nc = pixmap[np + 1];
	nc += quant_errorG * 1 / 16 * 255;
	pixmap[np + 1] = nc > 255 ? 255 : nc;
	nc = pixmap[np + 2];
	nc += quant_errorB * 1 / 16 * 255;
	pixmap[np + 2] = nc > 255 ? 255 : nc;

	return Color(newColorR, newColorG, newColorB);
}

Color burkesDither(int i, int j)
{
	// modification of floydSteinberg method
	// diffuses the error to 7 unprocessed neighboring pixels

	int pix = (i * width + j) * 3;
	int ny, nx, nc, np;
	float oldColorR, oldColorG, oldColorB;
	float newColorR, newColorG, newColorB;
	float quant_errorR, quant_errorG, quant_errorB;

	oldColorR = pixmap[pix] / 255.0;
	if (oldColorR < 0.25)
		newColorR = 0;
	else if (oldColorR >= 0.25 && oldColorR < 0.5)
		newColorR = 0.3;
	else if (oldColorR >= 0.5 && oldColorR < 0.75)
		newColorR = 0.7;
	else
		newColorR = 1.0;

	oldColorG = pixmap[pix + 1] / 255.0;
	if (oldColorG < 0.25)
		newColorG = 0;
	else if (oldColorG >= 0.25 && oldColorG < 0.5)
		newColorG = 0.3;
	else if (oldColorG >= 0.5 && oldColorG < 0.75)
		newColorG = 0.7;
	else
		newColorG = 1.0;

	oldColorB = pixmap[pix + 2] / 255.0;
	if (oldColorB < 0.25)
		newColorB = 0;
	else if (oldColorB >= 0.25 && oldColorB < 0.5)
		newColorB = 0.3;
	else if (oldColorB >= 0.5 && oldColorB < 0.75)
		newColorB = 0.7;
	else
		newColorB = 1.0;

	quant_errorR = oldColorR - newColorR;
	quant_errorG = oldColorG - newColorG;
	quant_errorB = oldColorB - newColorB;

	ny = i; nx = j + 1;
	np = (ny * width + nx) * 3;
	nc = pixmap[np];
	nc += quant_errorR * 8 / 32 * 255;
	pixmap[np] = nc > 255 ? 255 : nc;
	nc = pixmap[np + 1];
	nc += quant_errorG * 8 / 32 * 255;
	pixmap[np + 1] = nc > 255 ? 255 : nc;
	nc = pixmap[np + 2];
	nc += quant_errorB * 8 / 32 * 255;
	pixmap[np + 2] = nc > 255 ? 255 : nc;

	ny = i + 1; nx = j + 2;
	np = (ny * width + nx) * 3;
	nc = pixmap[np];
	nc += quant_errorR * 4 / 32 * 255;
	pixmap[np] = nc > 255 ? 255 : nc;
	nc = pixmap[np + 1];
	nc += quant_errorG * 4 / 32 * 255;
	pixmap[np + 1] = nc > 255 ? 255 : nc;
	nc = pixmap[np + 2];
	nc += quant_errorB * 4 / 32 * 255;
	pixmap[np + 2] = nc > 255 ? 255 : nc;

	ny = i + 1; nx = j - 2;
	np = (ny * width + nx) * 3;
	nc = pixmap[np];
	nc += quant_errorR * 2 / 32 * 255;
	pixmap[np] = nc > 255 ? 255 : nc;
	nc = pixmap[np + 1];
	nc += quant_errorG * 2 / 32 * 255;
	pixmap[np + 1] = nc > 255 ? 255 : nc;
	nc = pixmap[np + 2];
	nc += quant_errorB * 2 / 32 * 255;
	pixmap[np + 2] = nc > 255 ? 255 : nc;

	ny = i + 1; nx = j - 1;
	np = (ny * width + nx) * 3;
	nc = pixmap[np];
	nc += quant_errorR * 4 / 32 * 255;
	pixmap[np] = nc > 255 ? 255 : nc;
	nc = pixmap[np + 1];
	nc += quant_errorG * 4 / 32 * 255;
	pixmap[np + 1] = nc > 255 ? 255 : nc;
	nc = pixmap[np + 2];
	nc += quant_errorB * 4 / 32 * 255;
	pixmap[np + 2] = nc > 255 ? 255 : nc;

	ny = i + 1; nx = j;
	np = (ny * width + nx) * 3;
	nc = pixmap[np];
	nc += quant_errorR * 8 / 32 * 255;
	pixmap[np] = nc > 255 ? 255 : nc;
	nc = pixmap[np + 1];
	nc += quant_errorG * 8 / 32 * 255;
	pixmap[np + 1] = nc > 255 ? 255 : nc;
	nc = pixmap[np + 2];
	nc += quant_errorB * 8 / 32 * 255;
	pixmap[np + 2] = nc > 255 ? 255 : nc;

	ny = i + 1; nx = j + 1;
	np = (ny * width + nx) * 3;
	nc = pixmap[np];
	nc += quant_errorR * 4 / 32 * 255;
	pixmap[np] = nc > 255 ? 255 : nc;
	nc = pixmap[np + 1];
	nc += quant_errorG * 4 / 32 * 255;
	pixmap[np + 1] = nc > 255 ? 255 : nc;
	nc = pixmap[np + 2];
	nc += quant_errorB * 4 / 32 * 255;
	pixmap[np + 2] = nc > 255 ? 255 : nc;

	ny = i + 1; nx = j + 2;
	np = (ny * width + nx) * 3;
	nc = pixmap[np];
	nc += quant_errorR * 2 / 32 * 255;
	pixmap[np] = nc > 255 ? 255 : nc;
	nc = pixmap[np + 1];
	nc += quant_errorG * 2 / 32 * 255;
	pixmap[np + 1] = nc > 255 ? 255 : nc;
	nc = pixmap[np + 2];
	nc += quant_errorB * 2 / 32 * 255;
	pixmap[np + 2] = nc > 255 ? 255 : nc;

	return Color(newColorR, newColorG, newColorB);
}

Color randomDither(int i, int j)
{
	// randomly dither into 3 colors
	int np = (i * width + j) * 3;
	int nc = pixmap[np];
	nc += int(randomGen(0, 255));
	if (nc < 255 / 4)
		nc = 255 / 4;
	else if (nc < 255 / 2)
		nc = 255 / 2;
	else
		nc = 255;
	return Color(nc, nc, nc);
}

Color bayerDither(int i, int j)
{
	// Ordered dithering using Bayer matrices

	int np = (i * width + j) * 3;
	int nc = pixmap[np];
	 std::vector<std::vector<int>> D{ {3, 1}, {0, 2} };
	/*std::vector<std::vector<int>> D{	{ 1, 9, 3, 11 },
										{ 13, 5, 15, 7 },
										{ 4, 12, 2, 10 },
										{ 16, 8, 14, 6 }	};*/
	int n = D.size();
	int y = i % n;
	int x = j % n;
	int threshold = D[y][x];
	float r = 10;// n*n;
	
	float oldColorR, oldColorG, oldColorB;
	float newColorR, newColorG, newColorB;

	oldColorR = pixmap[np] / 255.0;
	oldColorR += float(threshold) / r;
	if (oldColorR < 0.25)
		newColorR = 0;
	else if (oldColorR >= 0.25 && oldColorR < 0.5)
		newColorR = 0.3;
	else if (oldColorR >= 0.5 && oldColorR < 0.75)
		newColorR = 0.7;
	else
		newColorR = 1.0;

	oldColorG = pixmap[np + 1] / 255.0;
	oldColorG += float(threshold) / r;
	if (oldColorG < 0.25)
		newColorG = 0;
	else if (oldColorG >= 0.25 && oldColorG < 0.5)
		newColorG = 0.3;
	else if (oldColorG >= 0.5 && oldColorG < 0.75)
		newColorG = 0.7;
	else
		newColorG = 1.0;

	oldColorB = pixmap[np + 2] / 255.0;
	oldColorB += float(threshold) / r;
	if (oldColorB < 0.25)
		newColorB = 0;
	else if (oldColorB >= 0.25 && oldColorB < 0.5)
		newColorB = 0.3;
	else if (oldColorB >= 0.5 && oldColorB < 0.75)
		newColorB = 0.7;
	else
		newColorB = 1.0;

	return Color(newColorR, newColorG, newColorB);
}

void render()
{
#ifdef COMPOSTION
	k.boxFilter();
#endif

#ifdef RASTERIZED_SHAPES
	width = 800;
	height = 800;
	pixmap.clear();
	pixmap.resize(width * height * 3);
#endif

#ifdef TRANSFORMATIONS
	// for bilinear warping
	std::vector<std::vector<double>> transFormedcorners;
	getPerspectiveTransformedCorners(transFormedcorners);
#endif

#pragma omp parellel for
	for (int i = 0; i < height; ++i)
	{
		for (int j = 0; j < width; ++j)
		{
			int p = (i * width + j) * 3;

#ifdef PROCEDURAL_IMAGE_GEN
			double y = double(i) / double(height);
			double x = double(j) / double(width);
			double r = randomGen();

			// 1.ppm
			pixmap[p + 0] = (1 + sin(10 * y)) * x * 0.5 * 255;
			pixmap[p + 1] = x * y * 255;
			pixmap[p + 2] = 0x60;

			// 2.ppm
			pixmap[p + 0] = x * 255;
			pixmap[p + 1] = (x * y) * 1.0 * 255;
			pixmap[p + 2] = (x + y) * 0.5 * 255;

			// 3.ppm
			if (i == 0 || j == 0)
				continue;
			int p_upLeft = ((i - 1) * width + (j - 1)) * 3;
			pixmap[p + 0] += pixmap[p_upLeft + 0] + r * x * 10 * 255;
			pixmap[p + 1] += pixmap[p_upLeft + 1] + r * y * 10 * 255;
			pixmap[p + 2] += pixmap[p_upLeft + 2] + r * x * y * 100 * 255;

			return;
#endif

#ifdef ANTIALIASING
			int subX = 3;
			int subY = 3;
			double rx = randomGen(0.0, 1.0 / subX);
			double ry = randomGen(0.0, 1.0 / subY);
#else
			int subX = 1;
			int subY = 1;
			double rx = 0.5;
			double ry = 0.5;
#endif
			double weighted = 1.0 / (subX * subY);
			double r = 0, g = 0, b = 0;
			double X, Y, x, y;

			for (int py = 0; py < subY; ++py)
			{
				for (int px = 0; px < subX; ++px)
				{
					Y = i + ry + (py + 0.5) / subY;
					X = j + rx + (px + 0.5) / subX;

					Color c;

#ifdef RASTERIZED_SHAPES
					c = shaded(X, Y);
					// c = shapeColor(star(X, Y));
#endif

#ifdef CONVOLUTION_FILTERS
					// Select any ONE filter type
					// Select any ONE filter within the type

					//ConvolutionFilter k(9, 9);
					//k.boxFilter();
					//k.radialFilter();
					//k.motionFilter();

					//MorphologicalFilter k(9, 9);
					//k.dilationFilter();
					//k.erosionFilter();

					DerivativeFilter k(3, 3);
					//k.embossFilter1();
					//k.embossFilter2();
					//k.embossFilter3();
					//k.embossFilter4();
					k.edgeFilter();

					c = k.eval(i, j);
#endif

#ifdef TRANSFORMATIONS
					//inverseTransform(X, Y);
					//bilinearWarpTransform(X, Y, transFormedcorners);
					warpTransform(X, Y);
#endif
					int pix = (int(Y) * width + int(X)) * 3;
#ifdef COMPOSITIONS
					if (pix < pixmap.size() && Y < height && X < width && Y >= 0 && X >= 0)
						c = compositeOperation(pix, i, j);
#endif

#ifdef DITHERING
					if (pix < pixmap.size() && int(Y) < height && int(X) < width && int(Y) >= 0 && int(X) >= 0)
					{
						// Enable any one of the following 

						// Random Dithering algorithm
						//c = randomDither(i, j);

						// each of the following uses 4 values per channel, that is, 64 color palette

						// Unordered or Error-Diffusion Dithering algorithms
						//c = floydSteinbergDither(i, j);
						//c = burkesDither(i, j);

						// Ordered Dithering algorithms
						c = bayerDither(i, j);
					}
#endif

					/*if (pix < pixmap.size() && int(Y) < height && int(X) < width && int(Y) >= 0 && int(X) >= 0)
						c = Color(pixmap[pix + 0] / 255.0, pixmap[pix + 1] / 255.0, pixmap[pix + 2] / 255.0);*/
					
					r += c.r * weighted;
					g += c.g * weighted;
					b += c.b * weighted;
				}
			}

			// store in a different pixmap; 
			// if you use same one, the next pixels will use convoluted pixel colors instead of input pixel colors
			pixmap2[p + 0] = (uint8_t)(r * 255);
			pixmap2[p + 1] = (uint8_t)(g * 255);
			pixmap2[p + 2] = (uint8_t)(b * 255);
		}
	}

	pixmap = std::move(pixmap2);
}