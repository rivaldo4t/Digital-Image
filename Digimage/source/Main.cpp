#include "main.hpp"

// First Part of Assignment
#define READ_WRITE_DISPLAY

// Second Part of Assignment
// #define PROCEDURAL_IMAGE_GEN

void render()
{
	width = 512;
	height = 512;
	pixmap.clear();
	pixmap.resize(width * height * 3);

	for (int i = 0; i < height; ++i)
	{
		for (int j = 0; j < width; ++j)
		{
			int p = (i * width + j) * 3;
			double y = double(i) / double(height);
			double x = double(j) / double(width);

			std::default_random_engine generator(428697);
			std::uniform_real_distribution<double> distribution(0.0, 1.0);
			double r = distribution(generator);

			// 1.ppm
			/*pixmap[p + 0] = (1 + sin(10*y)) * x * 0.5 * 255;
			pixmap[p + 1] = x * y * 255;
			pixmap[p + 2] = 0x60;*/

			// 2.ppm
			/*pixmap[p + 0] = x * 255;
			pixmap[p + 1] = (x * y) * 1.0 * 255;
			pixmap[p + 2] = (x + y) * 0.5 * 255;*/

			// 3.ppm
			if (i == 0 || j == 0)
				continue;
			int p_upLeft = ((i-1) * width + (j-1)) * 3;
			pixmap[p + 0] += pixmap[p_upLeft + 0] + r * x * 10 * 255;
			pixmap[p + 1] += pixmap[p_upLeft + 1] + r * y * 10 * 255;
			pixmap[p + 2] += pixmap[p_upLeft + 2] + r * x * y * 100 * 255;
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
	// pixmap = img.data; // displays flipped iamge with OpenGL
	int row = width * 3; // flip the data after reading
	for (int i = 0; i < img.data.size(); i += row)
		pixmap.insert(pixmap.begin(), img.data.begin() + i, img.data.begin() + i + row);
#endif

#ifdef PROCEDURAL_IMAGE_GEN
	render();
	Image img;
	img.w = width;
	img.h = height;
	img.data = pixmap;
	writePPM(img, "ppm/4.ppm");
	//writePseudoPPM(img, "ppm/a.txt");
#endif
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
	glutInitWindowPosition(600, 200);
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
