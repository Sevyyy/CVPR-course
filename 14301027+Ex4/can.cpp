#include <iostream>
#include <cmath>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include "can.h"
#include "CImg.h"
using namespace cimg_library;
using namespace std;


canny::canny(string file_path)
{
	// load the bmp and get the size
	img.load_bmp(file_path.c_str());
	width = img._width;
	height = img._height;

	//initial the several images
	gray_img = CImg<unsigned char>(width, height, 1, 1, 0);
	gaussian_img_x = CImg<float>(width, height, 1, 1, 0);
	gaussian_img_y = CImg<float>(width, height, 1, 1, 0);
	grad_x = CImg<float>(width, height, 1, 1, 0);
	grad_y = CImg<float>(width, height, 1, 1, 0);
	magnitude = CImg<int>(width, height, 1, 1, 0);
	edge = CImg<int>(width, height, 1, 1, 0);

	img.display();
}

void canny::to_gray_scale()
{
	cimg_forXY(img, x, y)
	{
		int r = img(x, y, 0);
		int g = img(x, y, 1);
		int b = img(x, y, 2);
		gray_img(x, y) = (unsigned char)(r * 0.2126 + g * 0.7152 + b * 0.0722);
	}
	gray_img.display();
	//gray_img.save(("../0.1/" + file_name + "_gray.bmp").c_str());
}

int canny::gaussian_filter(int kernel_width, int kernel_radius)
{
	//define the gaussian kernel and the diffrence kernel
	vector<float> kernel(kernel_width, 0);
	vector<float> diff_kernel(kernel_width, 0);

	// init the two kernel
	int kwidth;
	for(kwidth = 0; kwidth < kernel_width; kwidth++)
	{
		float g1, g2, g3;
		g1 = gaussian((float)kwidth, kernel_radius);
		if(g1 <= GAUSSIAN_CUT_OFF && kwidth >= 2)
			break;
		g2 = gaussian(kwidth-0.5f, kernel_radius);
		g3 = gaussian(kwidth+0.5f, kernel_radius);
		float temp = 2.0f * cimg::PI * kernel_radius * kernel_radius;
		kernel[kwidth] = (g1+g2+g3) / 3.0f / temp;
		diff_kernel[kwidth] = g3 - g2;
	}

	int init_x = kwidth - 1;
	int max_x = width - init_x;
	int init_y = init_x;
	int max_y = height - init_y;

	//do the converlution in x and y direction
	cimg_forXY(gray_img, x, y)
	{
		if(x >= init_x && x < max_x && y >= init_y && y < max_y)
		{
			float sum_x = gray_img(x,y) * kernel[0];
			float sum_y = sum_x;
			int w = 1;
			while(w < kwidth)
			{
				sum_y += kernel[w] * (gray_img(x,y-w) + gray_img(x,y+w));
				sum_x += kernel[w] * (gray_img(x-w,y) + gray_img(x+w,y));
				w++;
			}
			gaussian_img_x(x,y) = sum_x;
			gaussian_img_y(x,y) = sum_y;
		}
	}
	//gaussian_img_x.display();
	//gaussian_img_y.display();

	// get the gradient
	cimg_forXY(gray_img, x, y)
	{
		if(x >= init_x && x < max_x && y >= init_y && y < max_y)
		{
			float sum = 0.0f;
			for(int i = 0; i < kwidth; i++)
				sum += diff_kernel[i] * (gaussian_img_y(x-i,y) - gaussian_img_y(x+i,y));
			grad_x(x, y) = sum;
		}
	}

	cimg_forXY(gray_img, x, y)
	{
		if(x >= kwidth && x < width - kwidth && y >= init_y && y < max_y)
		{
			float sum = 0.0f;
			for(int i = 0; i < kwidth; i++)
				sum += diff_kernel[i] * (gaussian_img_y(x,y-i) - gaussian_img_y(x,y+i));
			grad_y(x, y) = sum;
		}
	}

	return kwidth;
}

void canny::non_maximal_supression(int kwidth)
{
	int init_x = kwidth;
	int max_x = width - kwidth;
	int init_y = kwidth;
	int max_y = height - kwidth;

	cimg_forXY(gray_img, x, y)
	{
		if(x >= init_x && x < max_x && y >= init_y && y < max_y)
		{
			float x_grad = grad_x(x,y);
			float y_grad = grad_y(x,y);
			float grad_mag = hypotenuse(x_grad, y_grad);

			//get the 8 direction magnitude
			float n_mag = hypotenuse(grad_x(x,y-1), grad_y(x,y-1));
			float s_mag = hypotenuse(grad_x(x,y+1), grad_y(x,y+1));
			float w_mag = hypotenuse(grad_x(x-1,y), grad_y(x-1,y));
			float e_mag = hypotenuse(grad_x(x+1,y), grad_y(x+1,y));
			float ne_mag = hypotenuse(grad_x(x+1,y-1), grad_y(x+1,y-1));
			float se_mag = hypotenuse(grad_x(x+1,y+1), grad_y(x+1,y+1));
			float sw_mag = hypotenuse(grad_x(x-1,y+1), grad_y(x-1,y+1));
			float nw_mag = hypotenuse(grad_x(x+1,y-1), grad_y(x+1,y-1));
			float tmp;

			//judge if is edge point or not(using interpolation)
			int flag = ((x_grad * y_grad <= 0.0f) 
					? ffabs(x_grad) >= ffabs(y_grad)
						? (tmp = ffabs(x_grad * grad_mag)) >= ffabs(y_grad * ne_mag - (x_grad + y_grad) * e_mag)
												   && tmp > ffabs(y_grad * sw_mag - (x_grad + y_grad) * w_mag)
						: (tmp = ffabs(y_grad * grad_mag)) >= ffabs(x_grad * ne_mag - (y_grad + x_grad) * n_mag)
												   && tmp > ffabs(x_grad * sw_mag - (y_grad + x_grad) * s_mag)
					: ffabs(x_grad) >= ffabs(y_grad)
						? (tmp = ffabs(x_grad * grad_mag)) >= ffabs(y_grad * se_mag + (x_grad - y_grad) * e_mag)
												   && tmp > ffabs(y_grad * nw_mag + (x_grad - y_grad) * w_mag)
						: (tmp = ffabs(y_grad * grad_mag)) >= ffabs(x_grad * se_mag + (y_grad - x_grad) * s_mag)
												   && tmp > ffabs(x_grad * nw_mag + (y_grad - x_grad) * n_mag)
					);
			if(flag)
				magnitude(x,y) = (grad_mag >= MAGNITUDE_LIMIT) ? MAGNITUDE_MAX : (int)(MAGNITUDE_SCALE * grad_mag);
			else
				magnitude(x,y) = 0;
		}
	}
	//magnitude.display();
	//magnitude.save(("../0.1/" + file_name + "_magnitude.bmp").c_str());
}

void canny::perform_hysteresis(int low, int high)
{
	int offset = 0;
	cimg_forXY(edge, x, y)
	{
		//process the still not detected strong point
		if(edge(x,y) == 0 && magnitude(x,y) >= high)
			follow(x, y, low);
	}
	edge.display();
	//edge.save(("../out/" + file_name + "_edge.bmp").c_str());
}

void canny::follow(int x1, int y1, int threshold)
{
	edge(x1, y1) = magnitude(x1, y1);
	int x0 = x1 == 0 ? x1 : x1-1;
	int x2 = x1 == width-1 ? x1 : x1+1;
	int y0 = y1 == 0 ? y1 : y1-1;
	int y2 = y1 == height-1 ? y1 : y1+1;

	for(int x = x0; x <= x2; x++)
	{
		for(int y = y0; y <= y2; y++)
		{
			//process the still not detected weak point
			if((y != y1 || x != x1) && edge(x,y) == 0 && magnitude(x,y) >= threshold)
				follow(x, y, threshold);
		}
	}
}

CImg<int> canny::get_edge_img()
{
	return edge;
}


float hypotenuse(float x, float y)
{
	return (float)sqrt(x*x + y*y);
}

float gaussian(float x, float sigma)
{
	return (float) exp(-(x*x) / (2.0f*sigma*sigma));
}