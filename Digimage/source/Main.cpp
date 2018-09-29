#include "main.hpp"
// #define READ_WRITE_DISPLAY
// #define PROCEDURAL_IMAGE_GEN

// to enable and diable anti aliasing
#define ANTIALIASING

// enable to see shading, disable to see other functions
#define SHADING

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
			pixmap[p + 0] = (1 + sin(10*y)) * x * 0.5 * 255;
			pixmap[p + 1] = x * y * 255;
			pixmap[p + 2] = 0x60;

			// 2.ppm
			pixmap[p + 0] = x * 255;
			pixmap[p + 1] = (x * y) * 1.0 * 255;
			pixmap[p + 2] = (x + y) * 0.5 * 255;

			// 3.ppm
			if (i == 0 || j == 0)
				continue;
			int p_upLeft = ((i-1) * width + (j-1)) * 3;
			pixmap[p + 0] += pixmap[p_upLeft + 0] + r * x * 10 * 255;
			pixmap[p + 1] += pixmap[p_upLeft + 1] + r * y * 10 * 255;
			pixmap[p + 2] += pixmap[p_upLeft + 2] + r * x * y * 100 * 255;
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
					
#ifndef SHADING
					// replace star() by function(), convex(), blobby() for other shapes
					if (star(X, Y) <= 0.0)
					{
						r += 226.0 / 256.0 * weighted;
						g += 48.0 / 256.0 * weighted;
						b += 81.0 / 256.0 * weighted;
					}
					else
					{
						r += 65.0 / 256.0 * weighted;
						g += 76.0 / 256.0 * weighted;
						b += 83.0 / 256.0 * weighted;
					}
#else
					Color c = shaded(X, Y);
					r += c.r * weighted;
					g += c.g * weighted;
					b += c.b * weighted;
#endif
				}
			}

			pixmap[p + 0] = (uint8_t)(r * 255);
			pixmap[p + 1] = (uint8_t)(g * 255);
			pixmap[p + 2] = (uint8_t)(b * 255);
		}
	}
}

void setPixels()
{

#ifdef READ_WRITE_DISPLAY
	// -------- Real PPM - Read / Write --------
	Image img = readPPM("ppm/construction.ppm");
	writePPM(img, "ppm/temp.ppm");

	// -------- Custom PPM - Read / Write --------
	// Image img = readPseudoPPM("ppm/pseudo.txt");
	// writePseudoPPM(img, "ppm/a.txt");

	// -------- Populate pixmap for display with OpenGL --------
	width = img.w;
	height = img.h;
	img.flip();
	pixmap = img.data;
#endif

	render();
	Image img(width, height, pixmap);
	img.flip();
	writePPM(img, "ppm/shaded.ppm");
}

static void windowResize(int w, int h)
{
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, (w / 2), 0, (h / 2), 0, 1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

static void windowDisplay(void)
{
	glClear(GL_COLOR_BUFFER_BIT);
	glRasterPos2i(0, 0);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glDrawPixels(width, height, GL_RGB, GL_UNSIGNED_BYTE, reinterpret_cast<char*>(pixmap.data()));
	glFlush();
}

static void processMouse(int button, int state, int x, int y)
{
	/*if (state == GLUT_UP)
		exit(0);*/
}

static void processKeyboard(unsigned char key, int x, int y)
{
	if (key == 27)
		exit(0);
}

static void init(void)
{
	glClearColor(1, 1, 1, 1);
}

int main(int argc, char *argv[])
{
	setPixels();

	glutInit(&argc, argv);
	glutInitWindowPosition(500, 50);
	glutInitWindowSize(width, height);
	glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
	glutCreateWindow("--------Digitized--------");
	init();
	glutReshapeFunc(windowResize);
	glutDisplayFunc(windowDisplay);
	glutMouseFunc(processMouse);
	glutKeyboardFunc(processKeyboard);
	glutMainLoop();

	return 0;
}
