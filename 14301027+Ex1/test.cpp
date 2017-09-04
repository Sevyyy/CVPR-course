#include <vector>
#include <algorithm>
#include <iostream>
#include "CImg.h"
using namespace cimg_library;
using namespace std;


void test_cimg()
{
	// read in the 1.bmp
	CImg<unsigned char> SrcImg;
	SrcImg.load_bmp("1.bmp");

	// display the image
	SrcImg.display();

	// define the two color
	const unsigned char blue[] = {0,0,255};
	const unsigned char yellow[] = {255,255,0};

	// make the transform
	cimg_forXY(SrcImg, x, y)
	{
		if (SrcImg(x, y, 0) == 255 && SrcImg(x, y, 1) == 255 && SrcImg(x, y, 2) == 255)
		{
			SrcImg(x, y, 1) = 0;
			SrcImg(x, y, 2) = 0;
		}
		if (SrcImg(x, y, 0) == 0 && SrcImg(x, y, 1) == 0 && SrcImg(x, y, 2) == 0)
		{
			SrcImg(x, y, 1) = 255;
		}
	}

	// draw the two circle
	SrcImg.draw_circle(50,50,30,blue);
	SrcImg.draw_circle(50,50,3,yellow);

	// display and save the image
	SrcImg.display();

	SrcImg.save("test.bmp");
}

int main(int argc, char* argv[])
{
	test_cimg();
	return 0;
}

