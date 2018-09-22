#include <iostream>
#include <fstream>
#include <GL/glut.h>
#include <main.hpp>
#include <memory>
#include <vector>
#include "image.hpp"
#include "ppm.hpp"

int width, height;
std::vector<uint8_t> pixmap;

void render()
{
	for (int i = 0; i < width; ++i)
	{
		for (int j = 0; j < height; ++j)
		{

		}
	}
}

void setPixels()
{
	// Real PPM
	// Image img = readPPM("ppm/cube.ppm");
	// writePPM(img, "ppm/cube2.ppm");

	// Custom PPM
	// Image img = readPseudoPPM("ppm/pseudo.txt");
	// writePseudoPPM(img, "ppm/a.txt");

	// width = img.w;
	// height = img.h;
	// pixmap = img.data;

	render();
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
	glDrawPixels(width, height, GL_RGB, GL_UNSIGNED_BYTE, pixmap.data());
	glFlush();
}

static void processMouse(int button, int state, int x, int y)
{
	if (state == GLUT_UP)
		exit(0);
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
	glutInitWindowPosition(400, 400);
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
