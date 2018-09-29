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
