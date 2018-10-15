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

// #define RASTERIZES_SHAPES
// #define PROCEDURAL_IMAGE_GEN
// #define ANTIALIASING

int width, height;
std::vector<uint8_t> pixmap;
std::vector<uint8_t> pixmap2;

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

	Color eval(int x, int y)
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

				int index = i * kWidth + j;

				c += neighborColor * weights[index];
			}
		}

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
		
		cGradY = Color(	(cGradY.r + 0.5*absWeightSum) / absWeightSum,
						(cGradY.g + 0.5*absWeightSum) / absWeightSum,
						(cGradY.b + 0.5*absWeightSum) / absWeightSum);

		cGradX = Color(	(cGradX.r + 0.5*absWeightSum2) / absWeightSum2,
						(cGradX.g + 0.5*absWeightSum2) / absWeightSum2,
						(cGradX.b + 0.5*absWeightSum2) / absWeightSum2);

		cGradMag = cGradY;
		cGradMag += cGradX;
		cGradMag = cGradMag * 0.5;
		cGradMag = Color(1 - cGradMag.r, 1 - cGradMag.g, 1 - cGradMag.b);

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

void render()
{
#ifdef RASTERIZES_SHAPES
	width = 800;
	height = 800;
	pixmap.clear();
	pixmap.resize(width * height * 3);
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
			int subX = 10;
			int subY = 10;
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
					y = Y / height;
					x = X / width;

#ifdef RASTERIZES_SHAPES
					Color c = shaded(X, Y);
#endif
					// Select any ONE filter type
					// Select any ONE filter within the type

					ConvolutionFilter k(9, 9);
					//k.boxFilter();
					//k.radialFilter();
					//k.motionFilter();

					//MorphologicalFilter k(9, 9);
					//k.dilationFilter();
					//k.erosionFilter();

					//DerivativeFilter k(3, 3);
					//k.embossFilter1();
					//k.embossFilter2();
					//k.embossFilter3();
					//k.embossFilter4();
					//k.edgeFilter();

					Color c = k.eval(i, j);

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