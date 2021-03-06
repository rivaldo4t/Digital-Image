#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include "image.hpp"

/*
	PseudoPPm format (ASCII txt file)
	----------------
	4T
	width height
	r0 b0 g0 r1 b1 c1 . . .(number of columns = (width * 3))
	.
	.
	.
	(number of rows = (height))
	----------------
*/

Image readPseudoPPM(std::string fileName)
{
	Image img;

	std::ifstream ifs(fileName.c_str(), std::ios::in);
	if (ifs.fail())
	{
		std::cerr << "Cannot open file " + fileName << std::endl;
		exit(1);
	}

	std::string header;
	ifs >> header;
	if (header != "4T")
	{
		std::cerr << "Header unsupported for file " + fileName << std::endl;
		exit(1);
	}

	int w, h;
	ifs >> w >> h;
	img.w = w;
	img.h = h;
	img.data.clear();
	img.data.reserve(w * h * 3);
	ifs.ignore(256, '\n');
	
	while (ifs.good())
	{
		int temp;
		ifs >> temp;
		img.data.push_back(temp);
	}

	ifs.close();
	std::cout << "Read successfully " << fileName << std::endl;
	return img;
}

void writePseudoPPM(Image img, std::string fileName)
{
	if (img.w < 0 || img.h < 0 || img.data.size() != img.w * img.h * 3)
	{
		std::cerr << "Invalid image data\n";
		exit(1);
	}

	std::ofstream ofs(fileName.c_str(), std::ios::out);
	ofs << "4T" << "\n"
		<< img.w << " " << img.h << "\n";

	for (unsigned int i = 0; i < img.h * img.w * 3; ++i)
	{
		ofs << (int)img.data[i];
		if ((i+1) % (img.w * 3) != 0)
			ofs << " ";
		else
			ofs << "\n";
	}

	ofs.close();
	std::cout << "Image written to " << fileName << std::endl;
}

Image readPPM(std::string fileName)
{
	Image img;

	std::ifstream ifs(fileName.c_str(), std::ios::binary);
	if (ifs.fail())
	{
		std::cerr << "Cannot open file " + fileName << std::endl;
		exit(1);
	}

	std::string header;
	ifs >> header;
	if (header != "P6")
	{
		std::cerr << "Header unsupported for file " + fileName << std::endl;
		exit(1);
	}

	// for skipping comments --------
	bool flag = true;
	char bufChar;
	char bufArray[100];
	ifs.get(bufChar); // skip '\n'
	while (flag)
	{
		ifs.get(bufChar); 
		if (bufChar != '#') // if this line is not a comment, put the char back in stream and continue reading file
		{
			flag = false;
			ifs.unget();
		}
		else
			ifs.getline(bufArray, 100, '\n'); // if this line is a comment, read the line and check for further comments in next iteration of loop
	}
	// ------------------------------

	int w, h, b;
	ifs >> w >> h >> b;
	img.w = w;
	img.h = h;
	img.data.resize(w * h * 3);
	ifs.ignore(256, '\n');
	ifs.read(reinterpret_cast<char *>(img.data.data()), img.data.size());

	ifs.close();
	std::cout << "Read successfully " << fileName << std::endl;
	return img;
}

void writePPM(Image img, std::string fileName)
{
	if (img.w < 0 || img.h < 0 || img.data.size() != img.w * img.h * 3)
	{
		std::cerr << "Invalid image data\n";
		exit(1);
	}

	std::ofstream ofs(fileName.c_str(), std::ios::binary);
	ofs << "P6" << "\n"
		<< img.w << " " << img.h << "\n"
		<< "255" << "\n";

	ofs.write(reinterpret_cast<char const*>(img.data.data()), img.data.size());
	ofs.close();
	std::cout << "Image written to " << fileName << std::endl;
}
