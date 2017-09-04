#include <iostream>
#include <cmath>
#include <vector>
#include "CImg.h"

using namespace std;
using namespace cimg_library;

// some math function
double gamma(double x);
double f(double x);
double ff(double x);
double mean(CImg<double> src);
double standard_deviations(CImg<double> src);

void to_gray_scale(CImg<unsigned char> &src_img, CImg<unsigned char> &gray_img);
vector<int> get_histogram(CImg<unsigned char> img);
CImg<unsigned char> histogram_equalization(CImg<unsigned char> img);
CImg<unsigned char> color_histogram_equalization(CImg<unsigned char> src);
CImg<double> rgb2xyz(CImg<double> src);
CImg<double> xyz2lab(CImg<double> src);
CImg<double> rgb2lab(CImg<double> src);
CImg<double> lab2xyz(CImg<double> src);
CImg<double> xyz2rgb(CImg<double> src);
CImg<double> lab2rgb(CImg<double> src);
CImg<double> color_transfer(CImg<double> src, CImg<double> tar);

int main()
{
	string file_path1 = "./image/1/";
	string file_path2 = "./image/2/";

	//1.histogram equalization test(including gray image and colorful image)
	for(int i = 0;i < 10;i++)
	{
		string file = file_path1 + "000" + (char)(i+'0') + ".bmp";
		CImg<unsigned char> src;
		src.load_bmp(file.c_str());
		CImg<unsigned char> gray(src._width, src._height, 1, 1, 0);
		to_gray_scale(src, gray);
		CImg<unsigned char> gray_he = histogram_equalization(gray);
		CImg<unsigned char> src_he = color_histogram_equalization(src);
		gray.display();
		gray_he.display();
		src.display();
		src_he.display();
	}

	//2.color transfer test(of course for color image)
	for(int i = 0;i < 10;i++)
	{
		string src_file = file_path2 + "100" + (char)(i+'0') + ".bmp";
		string tar_file = file_path2 + "200" + (char)(i+'0') + ".bmp";
		CImg<unsigned char> src, tar;
		src.load_bmp(src_file.c_str());
		tar.load_bmp(tar_file.c_str());
		CImg<double> out = color_transfer(src, tar);
		src.display();
		tar.display();
		out.display();
	}
	
	return 0;
}

double gamma(double x)
{
	return x > 0.04045 ? pow((x+0.055)/1.055,2.4) : x/12.92;
}

double f(double x)
{
	return x > 0.008856 ? pow(x, 1.0/3.0) : 7.787*x+16.0/116;
}

double ff(double x)
{
	return x > 0.206893? x*x*x : (x-16.0/116)/7.787;
}

void to_gray_scale(CImg<unsigned char> &src_img, CImg<unsigned char> &gray_img)
{
	cimg_forXY(src_img, x, y)
	{
		int r = src_img(x, y, 0);
		int g = src_img(x, y, 1);
		int b = src_img(x, y, 2);
		gray_img(x, y) = (unsigned char)(r*0.2126 + g*0.7152 + b*0.0722);
	}
}

vector<int> get_histogram(CImg<unsigned char> img)
{
	vector<int> h(256,0);
	cimg_forXY(img, x, y)
	{
		h[img(x,y)]++;
	}
	return h;
}

CImg<unsigned char> histogram_equalization(CImg<unsigned char> img)
{
	int width = img._width;
	int height = img._height;
	CImg<unsigned char> out(width, height, 1, 1, 0);
	vector<int> h = get_histogram(img);
	vector<int> p;
	int s = 0;
	for(int i = 0;i < 256; i++)
	{
		s += h[i];
		p.push_back((int)(256*s/(width*height)));
	}

	//perform the equalization
	cimg_forXY(img, x, y)
	{
		out(x,y) = p[img(x,y)] > 255 ? 255 : p[img(x,y)];
	}
	return out;
}

CImg<unsigned char> color_histogram_equalization(CImg<unsigned char> src)
{
	CImg<unsigned char> out(src._width, src._height, 1, 3, 0);

	CImg<unsigned char> r = histogram_equalization(src.get_channel(0));
	CImg<unsigned char> g = histogram_equalization(src.get_channel(1));
	CImg<unsigned char> b = histogram_equalization(src.get_channel(2));

	// operate the 3 channel respectively
	cimg_forXY(out, x, y)
	{
		out(x, y, 0) = r(x, y);
		out(x, y, 1) = g(x, y);
		out(x, y, 2) = b(x, y);
	}
	return out;
}

CImg<double> rgb2xyz(CImg<double> src)
{
	CImg<double> out(src._width, src._height, 1, 3, 0);
	cimg_forXY(src, x, y)
	{
		double R = src(x,y,0)*1/255;
		double G = src(x,y,1)*1/255;
		double B = src(x,y,2)*1/255;
		out(x,y,0) = R*0.414253 + G*0.357580 + B*0.180423;
		out(x,y,1) = R*0.212671 + G*0.715160 + B*0.072169;
		out(x,y,2) = R*0.019334 + G*0.119193 + B*0.950227;
	}
	return out;
}

CImg<double> xyz2lab(CImg<double> src)
{
	CImg<double> out(src._width, src._height, 1, 3, 0);
	double xn = 0.414253+0.357580+0.180423;
	double yn = 0.212671+0.715160+0.072169;
	double zn = 0.019334+0.119196+0.950227;
	cimg_forXY(src, x, y)
	{	
		double X = src(x,y,0) / xn;
		double Y = src(x,y,1) / yn;
		double Z = src(x,y,2) / zn;
		double fx = f(X);
		double fy = f(Y);
		double fz = f(Z);
		out(x,y,0) = cimg::max(0, 116*fy-16);
		out(x,y,1) = 500.0 * (fx - fy);
		out(x,y,2) = 200.0 * (fy - fz);
	}
	return out;
}

CImg<double> rgb2lab(CImg<double> src)
{
	return xyz2lab(rgb2xyz(src));
}

CImg<double> lab2xyz(CImg<double> src)
{
	CImg<double> out(src._width, src._height, 1, 3, 0);
	double Xn = 0.412453 + 0.357580 + 0.180423;
    double Yn = 0.212671 + 0.715160 + 0.072169;
    double Zn = 0.019334 + 0.119193 + 0.950227;
    cimg_forXY(src, x, y)
    {
    	double L = src(x,y,0);
    	double a = src(x,y,1);
    	double b = src(x,y,2);
    	double cY = (L+16)/116.0;
    	out(x,y,1) = Yn * ff(cY);
    	double cX = a/500 + cY;
    	out(x,y,0) = Xn * ff(cX);
    	double cZ = cY - b/200;
    	out(x,y,2) = Zn * ff(cZ);
    }
    return out;
}

CImg<double> xyz2rgb(CImg<double> src)
{
	CImg<double> out(src._width, src._height, 1, 3, 0);
	double Xn = 0.412453 + 0.357580 + 0.180423;
    double Yn = 0.212671 + 0.715160 + 0.072169;
    double Zn = 0.019334 + 0.119193 + 0.950227;
    cimg_forXY(src, x, y)
    {
    	double X = src(x,y,0) * 255;
    	double Y = src(x,y,1) * 255;
    	double Z = src(x,y,2) * 255;
    	double R = 3.240479*X  - 1.537150*Y - 0.498535*Z;
        double G = -0.969256*X + 1.875992*Y + 0.041556*Z;
        double B = 0.055648*X  - 0.204043*Y + 1.057311*Z;
        out(x,y,0) = R < 0 ? 0 : R > 255 ? 255 : R;
        out(x,y,1) = G < 0 ? 0 : G > 255 ? 255 : G;
        out(x,y,2) = B < 0 ? 0 : B > 255 ? 255 : B;
    }
    return out;
}

CImg<double> lab2rgb(CImg<double> src)
{
	return xyz2rgb(lab2xyz(src));
}

double mean(CImg<double> src)
{
	double sum = 0;
	cimg_forXY(src, x, y)
	{
		sum += src(x,y);
	}
	return sum / src.size();
}

double standard_deviations(CImg<double> src)
{
	double m = mean(src);
	double sum = 0;
	cimg_forXY(src, x, y)
	{
		sum += (src(x,y)-m)*(src(x,y)-m);
	}
	return sqrt(sum / src.size());
}

CImg<double> color_transfer(CImg<double> src, CImg<double> tar)
{
	CImg<double> out(src._width, src._height, 1, 3, 0);

	//transform to lab space
	CImg<double> src_lab = rgb2lab(src);
	CImg<double> tar_lab = rgb2lab(tar);

	// get the mean of all channel
	double sm_l = mean(src_lab.get_channel(0));
	double sm_a = mean(src_lab.get_channel(1));
	double sm_b = mean(src_lab.get_channel(2));
	double tm_l = mean(tar_lab.get_channel(0));
	double tm_a = mean(tar_lab.get_channel(1));
	double tm_b = mean(tar_lab.get_channel(2));

	//calculate the parameter
	double para_l = standard_deviations(tar_lab.get_channel(0))
			/standard_deviations(src_lab.get_channel(0));
	double para_a = standard_deviations(tar_lab.get_channel(1))
			/standard_deviations(src_lab.get_channel(1));
	double para_b = standard_deviations(tar_lab.get_channel(2))
			/standard_deviations(src_lab.get_channel(2));

	//perform the transform
	cimg_forXY(out,x,y)
	{
		out(x,y,0) = para_l * (src_lab(x,y,0) - sm_l) + tm_l;
		out(x,y,1) = para_a * (src_lab(x,y,1) - sm_a) + tm_a;
		out(x,y,2) = para_b * (src_lab(x,y,2) - sm_b) + tm_b;
	}

	//back to rgb space
	out = lab2rgb(out);
	return out;
}