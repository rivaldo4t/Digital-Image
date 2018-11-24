#include "main.hpp"

void setPixels()
{
	// make sure all images have same size
	Image img = readPPM("ppm/carv2.ppm");
	width = img.w;
	height = img.h;
	img.flip();
	pixmap = img.data;
	pixmap2 = pixmap; // resize instead of copy
	
	/*Image img2 = readPPM("ppm/dp_green.ppm");
	img2.flip();
	pixmap2 = img2.data;

	Image img3 = readPPM("ppm/dp_alpha.ppm");
	img3.flip();
	pixmap3 = img3.data;*/

	//render();

	/*Image outputImg(width, height, pixmap);
	outputImg.flip();
	writePPM(outputImg, "ppm/dither_o.ppm");*/
	

	// seam carving
	int deleteNumPixels = 200;
	for (int w = 0; w < deleteNumPixels; w++)
	{
		calculateEnergy(w);
		calculateMinEnergy(w);
		deleteSeam(w);
		
		Image outputImg1(width, height, pixmap2);
		outputImg1.flip();
		writePPM(outputImg1, "ppm/seam_" + std::to_string(2*w) + ".ppm");

		Image outputImg2(width, height, pixmap);
		outputImg2.flip();
		writePPM(outputImg2, "ppm/seam_"+ std::to_string(2*w+1) + ".ppm");
	}
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
