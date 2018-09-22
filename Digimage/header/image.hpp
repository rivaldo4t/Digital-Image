#pragma once
#include <iostream>
#include <vector>

class Image
{
public:
	unsigned int w, h;
	std::vector<uint8_t> data;
	Image() : w(0), h(0), data(std::vector<uint8_t>()) {}
};
