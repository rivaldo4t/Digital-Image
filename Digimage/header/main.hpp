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

int width, height;
std::vector<uint8_t> pixmap;