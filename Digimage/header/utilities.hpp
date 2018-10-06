#pragma once
#include <iostream>
#include <fstream>
#include <GL/glut.h>
#include <memory>
#include <vector>
#include <numeric>
#include <random>
#include <chrono>
#include "image.hpp"
#include "ppm.hpp"
#define _USE_MATH_DEFINES
#include <math.h>

// #define PROCEDURAL_IMAGE_GEN
// #define ANTIALIASING

int width, height;
std::vector<uint8_t> pixmap;

struct Color
{
	double r = 0, g = 0, b = 0;
};

template <class T>
T clamp(T val, T min = 0.0, T max = 1.0)
{
	if (val < min)
		val = min;
	else if (val > max)
		val = max;
	return val;
}

void RGBtoHSV(int r, int g, int b, double &h, double &s, double &v) 
{

#define maximum(x, y, z) ((x) > (y)? ((x) > (z)? (x) : (z)) : ((y) > (z)? (y) : (z))) 
#define minimum(x, y, z) ((x) < (y)? ((x) < (z)? (x) : (z)) : ((y) < (z)? (y) : (z)))

	double red, green, blue;
	double max, min, delta;

	red = r / 255.0; green = g / 255.0; blue = b / 255.0;  /* r, g, b to 0 - 1 scale */

	max = maximum(red, green, blue);
	min = minimum(red, green, blue);

	v = max;        /* value is maximum of r, g, b */

	if (max == 0) {    /* saturation and hue 0 if value is 0 */
		s = 0;
		h = 0;
	}
	else {
		s = (max - min) / max;           /* saturation is color purity on scale 0 - 1 */

		delta = max - min;
		if (delta == 0)                    /* hue doesn't matter if saturation is 0 */
			h = 0;
		else {
			if (red == max)                  /* otherwise, determine hue on scale 0 - 360 */
				h = (green - blue) / delta;
			else if (green == max)
				h = 2.0 + (blue - red) / delta;
			else /* (blue == max) */
				h = 4.0 + (red - green) / delta;
			h = h * 60.0;
			if (h < 0)
				h = h + 360.0;
		}
	}
}

void HSVtoRGB(int& R, int& G, int& B, double h, double s, double v) 
{
	double hh, p, q, t, ff;
	long i;
	double r, g, b;
	if (s <= 0.0) // < is bogus, just shuts up warnings
	{
		r = v;
		g = v;
		b = v;
		R = r * 255;
		G = g * 255;
		B = b * 255;
		return;
	}

	hh = h;
	if (hh >= 360.0) hh = 0.0;
	hh /= 60.0;
	i = (long)hh;
	ff = hh - i;
	p = v * (1.0 - s);
	q = v * (1.0 - (s * ff));
	t = v * (1.0 - (s * (1.0 - ff)));

	switch (i) 
	{
		case 0:
			r = v;
			g = t;
			b = p;
			break;
		case 1:
			r = q;
			g = v;
			b = p;
			break;
		case 2:
			r = p;
			g = v;
			b = t;
			break;
		case 3:
			r = p;
			g = q;
			b = v;
			break;
		case 4:
			r = t;
			g = p;
			b = v;
			break;
		case 5:
		default:
			r = v;
			g = p;
			b = q;
			break;
	}
	R = r * 255;
	G = g * 255;
	B = b * 255;
}

double randomGen(double min = 0.0, double max = 1.0)
{
	std::default_random_engine generator(428697);
	std::uniform_real_distribution<double> distribution(min, max);
	return distribution(generator);
}

double convex(double x, double y)
{
	double f = 2 * (x - 100) - y;
	double g = (x + 200) - 2 * y;
	double h = 400 - y;
	double k = 400 - x;
	return std::max(std::max(-f, g), std::max(-h, -k));
}

double function(double x, double y)
{
	return -(100 * sin(x / 50) + x - y);
}

double blobby(double x, double y)
{
	double f = 100 * sin(x / 50 + 400) + x - y;
	double g = 200 * cos(x / 100 - 70) + 250 - y;
	double h = 100 * sin(x / 200) - y;
	double k = 400 * cos(x / 400) - y;
	return std::max(std::max(-f, g), std::max(h, -k));
}

double star(double x, double y)
{
	double f = x + 4000 - 10 * y;
	double g = -4 * x + 4500 - 10 * y;
	double h = -2 * x + 1400 - y;
	double k = 2 * x - 200 - y;
	double m = 5 * x + 1000 - 10 * y;
	double t1 = std::max(std::max(-f, g), m);
	double t2 = std::max(std::max(-h, -k), m);
	double t3 = std::max(std::max(g, -h), -k);
	return std::min(std::min(t1, t2), t3);
}

Color shapeColor(double functionVal)
{
	Color c;
	if (functionVal <= 0.0)
	{
		c.r += 226.0 / 256.0;
		c.g += 48.0 / 256.0;
		c.b += 81.0 / 256.0;
	}	
	else
	{
		c.r += 65.0 / 256.0;
		c.g += 76.0 / 256.0;
		c.b += 83.0 / 256.0;
	}
	return c;
}

Color shaded(double x, double y)
{
	Color c;
	double f = pow(x - 400, 2) + pow(y - 400, 2) - pow(200, 2);
	double g = pow(x - 350, 2) + pow(y - 350, 2) - pow(100, 2);
	if (f <= 0)
	{
		c.r = (0 - f / pow(200, 2)) * 226.0 / 256.0;
		c.g = (0 - f / pow(200, 2)) * 48.0 / 256.0;
		c.b = 81.0 / 256.0;
		if (g <= 0)
		{
			c.r += pow((0 - g / pow(100, 2)), 2);
			c.g += pow((0 - g / pow(100, 2)), 2);
			c.b += pow((0 - g / pow(100, 2)), 2);
		}
	}
	else
	{
		c.r = 65.0 / 256.0;
		c.g = 76.0 / 256.0;
		c.b = 83.0 / 256.0;
	}
	c.r = clamp<double>(c.r);
	c.g = clamp<double>(c.g);
	c.b = clamp<double>(c.b);
	
	return c;
}

void render()
{
	width = 800;
	height = 800;
	pixmap.clear();
	pixmap.resize(width * height * 3);

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

					// replace star() by function(), convex(), blobby() or shaded() for other shapes
					//Color c = shaded(X, Y);
					Color c = shapeColor(star(X, Y));
					r += c.r * weighted;
					g += c.g * weighted;
					b += c.b * weighted;
				}
			}

			pixmap[p + 0] = (uint8_t)(r * 255);
			pixmap[p + 1] = (uint8_t)(g * 255);
			pixmap[p + 2] = (uint8_t)(b * 255);
		}
	}
}

// stored as (r0, r1, r2, r3 ..., rn)
// points are given as (0, r0), (1/n, r1), (2/n, r2), (3/n, r3) and (1, rn)
double pieceWiseLinearInterpolation(std::vector<double>& curveParam, double c)
{
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