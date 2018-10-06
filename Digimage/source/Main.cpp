#include "utilities.hpp"
// #define READ_WRITE_DISPLAY

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

	//render();

	
	// Project 3.1 - Curves Adjustment using piecewise linear interpolation
	// manipulate the img and store result in pixmap
#if 0
	 Image img = readPPM("ppm/ca.ppm");
	 curvesAdjustment(img);
	
	 // write the manipulated pixmap to file
	 Image img1(width, height, pixmap);
	 img1.flip();
	 writePPM(img1, "ppm/ca_m.ppm");
#endif

	// Project 3.2 - Hue replacement
	// Replace hue in destination image from hue in source image
#if 1
	Image src = readPPM("ppm/hue_src_7.ppm");
	Image dst = readPPM("ppm/hue_dst.ppm");
	replaceHues(src, dst);

	// write the manipulated pixmap to file
	Image img2(width, height, pixmap);
	img2.flip();
	writePPM(img2, "ppm/hue_replace_3.ppm");
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
