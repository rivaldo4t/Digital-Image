#pragma once
#include <iostream>
#include <vector>

class Image
{
public:
	unsigned int w, h;
	std::vector<uint8_t> data;
	Image() : w(0), h(0), data(std::vector<uint8_t>()) {}
	Image(unsigned int width, unsigned int height, std::vector<uint8_t> imgData = std::vector<uint8_t>()) : 
		w(width), h(height), data(imgData) {}

	// flip the data upside down
	void flip()
	{
		std::vector<uint8_t> flippedData;
		int row = w * 3;
		for (size_t i = 0; i < data.size(); i += row)
			flippedData.insert(flippedData.begin(), data.begin() + i, data.begin() + i + row);
		data = std::move(flippedData);
	}
};
