#pragma once
#include <iostream>
#include <fstream>
#include <GL/glut.h>
#include <main.hpp>
#include <memory>
#include <vector>
#include <numeric>
#include <random>
#include <chrono>
#include "image.hpp"
#include "ppm.hpp"
#define _USE_MATH_DEFINES
#include <math.h>

int width, height;
std::vector<uint8_t> pixmap;

struct Color
{
	double r = 0, g = 0, b = 0;
};

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
	if (c.r > 1)
		c.r = 1;
	if (c.g > 1)
		c.g = 1;
	if (c.b > 1)
		c.b = 1;
	return c;
}