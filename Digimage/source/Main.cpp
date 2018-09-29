#include "main.hpp"
// #define READ_WRITE_DISPLAY
// #define PROCEDURAL_IMAGE_GEN

// to enable and diable anti aliasing
#define ANTIALIASING

// enable to see shading, disable to see other functions
#define SHADING

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
					// replace star() by function(), convex() or blobby() for other shapes
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
	writePPM(img, "ppm/x.ppm");
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
