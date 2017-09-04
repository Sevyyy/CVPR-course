#ifndef _CANNY_H_
#define _CANNY_H_

#include "CImg.h"
#include <cstring>
#include <string>

using namespace cimg_library;
using namespace std;

#define ffabs(x) ( (x) >= 0 ? (x) : -(x) )
#define GAUSSIAN_CUT_OFF 0.005f
#define MAGNITUDE_SCALE 100.0f
#define MAGNITUDE_LIMIT 1000.0f
#define MAGNITUDE_MAX ((int) (MAGNITUDE_SCALE * MAGNITUDE_LIMIT))

class canny
{
private:
	string file_name;
	CImg<unsigned char> img;          //the origin image
	int width;
	int height;
	CImg<unsigned char> gray_img;     //the gray scale image
	CImg<float> gaussian_img_x;       //filtered by gaussian in x direction
	CImg<float> gaussian_img_y;       //filtered by gaussian in y direction
	CImg<float> grad_x;               //gradient in x direction
	CImg<float> grad_y;               //gradient in y direction
	CImg<int> magnitude;              //the edge magnitude
	CImg<int> edge;                   //the edge detected
public:
	canny(string file_path, string file_name);                                        //constructor
	void to_gray_scale();                                      //rgb to gray scale
	int gaussian_filter(int kernel_width, int kernel_radius);  //do the gaussian filter
	void non_maximal_supression(int kwidth);                   //do the non-maximal supression
	void perform_hysteresis(int low, int high);                //perform hysteresis
	void follow(int x1, int y1, int threshold);                //the helper in hysteresis
};

float hypotenuse(float x, float y);

float gaussian(float x, float sigma);



#endif