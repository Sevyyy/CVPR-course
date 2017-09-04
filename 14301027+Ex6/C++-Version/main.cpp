#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <set>
#include <algorithm>
#include "CImg.h"
#include "Sift.h"

using namespace std;
using namespace cimg_library;

#define PI 3.1415926
#define ANGLE 15.0
#define MATCH_THRESH 0.7

struct match_pair
{
	int i;  //index of img1
	int j;  //index of img2
	double dis;  //distance of two discriptor
};

// do interpolation
template <class T>
T interpolation(const CImg<T>& image, float x, float y, int channel); 
//do the cylinder projection
CImg<double> get_cylinder(const CImg<double> &src);
//get gray scale
CImg<double> get_gray_scale(CImg<double> img);
//get distance of two vector
double get_distance(vector<float> v1, vector<float> v2);
//get H
vector<double> get_transform_matrix(vector<pair<double, double> > uv, vector<pair<double,double> > xy);
//transform point base on H
pair<double, double> point_transform(vector<double> para, int u, int v);
//get inverse matrix
vector<double> inverse(vector<double> matrix);
//do sift-match
vector<match_pair> sift_match(const vector<SiftDescriptor> &desc1, const vector<SiftDescriptor> &desc2);
//do ransac
vector<double> ransac(const vector<match_pair> &matches, const vector<SiftDescriptor> &desc1, const vector<SiftDescriptor> &desc2);
//do stitch
CImg<double> stitch(vector<double> trans, const CImg<double> &img1, const CImg<double> &img2);
//whole image stitch
CImg<double> img_stitch(const CImg<double> &img1, const CImg<double> &img2);
//judge if a pixel is black
bool isEmpty(const CImg<unsigned char> &img, int x, int y);
//blend two image using pyramid blending algorithm
CImg<double> blendTwoImages(const CImg<double> &a, const CImg<double> &b);

int main()
{	
	CImg<double> cur("./input_s/C00.JPG");
	cur = get_cylinder(cur);
	cur.save("cur.jpg");
	CImg<double> next("./input_s/C01.JPG");
	next = get_cylinder(next);
	cur = img_stitch(cur,next);
	cur.save("cur.jpg");
	next = CImg<double>("./input_s/C02.JPG");
	next = get_cylinder(next);
	cur = img_stitch(cur,next);
	cur.save("cur.jpg");
	next = CImg<double>("./input_s/C03.JPG");
	next = get_cylinder(next);
	cur = img_stitch(cur,next);
	cur.save("cur.jpg");
	next = CImg<double>("./input_s/C04.JPG");
	next = get_cylinder(next);
	cur = img_stitch(cur,next);
	cur.save("cur.jpg");
	/*
	next = CImg<double>("./input_s/C05.JPG");
	next = get_cylinder(next);
	cur = img_stitch(cur,next);
	cur.save("cur.jpg");
	next = CImg<double>("./input_s/C06.JPG");
	next = get_cylinder(next);
	cur = img_stitch(cur,next);
	cur.save("cur.jpg");
	next = CImg<double>("./input_s/C07.JPG");
	next = get_cylinder(next);
	cur = img_stitch(cur,next);
	cur.save("cur.jpg");
	next = CImg<double>("./input_s/C08.JPG");
	next = get_cylinder(next);
	cur = img_stitch(cur,next);
	cur.save("cur.jpg");
	next = CImg<double>("./input_s/C09.JPG");
	next = get_cylinder(next);
	cur = img_stitch(cur,next);
	cur.save("cur.jpg");
	next = CImg<double>("./input_s/C10.JPG");
	next = get_cylinder(next);
	cur = img_stitch(cur,next);
	cur.save("cur.jpg");
	next = CImg<double>("./input_s/C11.JPG");
	next = get_cylinder(next);
	cur = img_stitch(cur,next);
	cur.save("cur.jpg");
	next = CImg<double>("./input_s/C12.JPG");
	next = get_cylinder(next);
	cur = img_stitch(cur,next);
	cur.save("cur.jpg");
	next = CImg<double>("./input_s/C13.JPG");
	next = get_cylinder(next);
	cur = img_stitch(cur,next);
	cur.save("cur.jpg");
	next = CImg<double>("./input_s/C14.JPG");
	next = get_cylinder(next);
	cur = img_stitch(cur,next);
	cur.save("cur.jpg");
	next = CImg<double>("./input_s/C15.JPG");
	next = get_cylinder(next);
	cur = img_stitch(cur,next);
	cur.save("cur.jpg");
	next = CImg<double>("./input_s/C16.JPG");
	next = get_cylinder(next);
	cur = img_stitch(cur,next);
	cur.save("cur.jpg");
	next = CImg<double>("./input_s/C17.JPG");
	next = get_cylinder(next);
	cur = img_stitch(cur,next);
	cur.save("cur.jpg");
	*/
	return 0;
}

template <class T>
T interpolation(const CImg<T>& image, float x, float y, int channel) 
{
    int x_pos = floor(x);
    float x_u = x - x_pos;
    int xb = (x_pos < image.width() - 1) ? x_pos + 1 : x_pos;

    int y_pos = floor(y);
	float y_v = y - y_pos;
    int yb = (y_pos < image.height() - 1) ? y_pos + 1 : y_pos;

	float P1 = image(x_pos, y_pos, channel) * (1 - x_u) + image(xb, y_pos, channel) * x_u;
	float P2 = image(x_pos, yb, channel) * (1 - x_u) + image(xb, yb, channel) * x_u;

    return P1 * (1 - y_v) + P2 * y_v;
}

//0.9840 >= ans >= 0.0160
CImg<double> get_cylinder(const CImg<double> &src)
{
	int projection_width, projection_height;
	CImg<double> res(src.width(), src.height(), 1, src.spectrum(), 0);
	float r;
	projection_width = src.width();
	projection_height = src.height();

	r = (projection_width / 2.0) / tan(ANGLE * PI / 180.0);

	for (int i = 0; i < res.width(); i++) 
	{
		for (int j = 0; j < res.height(); j++) 
		{
			float dst_x = i - projection_width / 2;
			float dst_y = j - projection_height / 2;

			float k = r / sqrt(r * r + dst_x * dst_x);
			float src_x = dst_x / k;
			float src_y = dst_y / k;

			if (src_x + projection_width / 2 >= 0 && src_x + projection_width / 2 < src.width()
				&& src_y + projection_height / 2 >= 0 && src_y + projection_height / 2 < src.height()) {
				for (int k = 0; k < res.spectrum(); k++) {
					res(i, j, k) = interpolation(src, src_x + projection_width / 2, src_y + projection_height / 2, k);
				}
			}
		}
	}
	return res;
}

CImg<double> get_gray_scale(CImg<double> img)
{
	CImg<double> gray(img.width(), img.height(), 1, 1, 0);
	cimg_forXY(img, x, y)
	{
		int r = img(x, y, 0);
		int g = img(x, y, 1);
		int b = img(x, y, 2);
		gray(x, y) = (r * 0.2126 + g * 0.7152 + b * 0.0722);
	}
	return gray;
}

double get_distance(vector<float> v1, vector<float> v2)
{
	int n = v1.size();
	double ans = 0;
	for(int i = 0;i < n;i++)
	{
		ans += (v1[i]-v2[i])*(v1[i]-v2[i]);
	}
	ans = sqrt(ans);
	return ans;
}

vector<double> get_transform_matrix(vector<pair<double, double> > uv, vector<pair<double,double> > xy)
{
	//get the 8 point
	double u1 = uv[0].first;
	double u2 = uv[1].first;
	double u3 = uv[2].first;
	double u4 = uv[3].first;
	double x1 = xy[0].first;
	double x2 = xy[1].first;
	double x3 = xy[2].first;
	double x4 = xy[3].first;
	double v1 = uv[0].second;
	double v2 = uv[1].second;
	double v3 = uv[2].second;
	double v4 = uv[3].second;
	double y1 = xy[0].second;
	double y2 = xy[1].second;
	double y3 = xy[2].second;
	double y4 = xy[3].second;

	//init the A from the 8 point
	double A[8][9] = {
					{x1, y1, 1, 0, 0, 0, -u1*x1, -u1*y1, u1},
					{0, 0, 0, x1, y1, 1, -v1*x1, -v1*y1, v1},
					{x2, y2, 1, 0, 0, 0, -u2*x2, -u2*y2, u2},
					{0, 0, 0, x2, y2, 1, -v2*x2, -v2*y2, v2},
					{x3, y3, 1, 0, 0, 0, -u3*x3, -u3*y3, u3},
					{0, 0, 0, x3, y3, 1, -v3*x3, -v3*y3, v3},
					{x4, y4, 1, 0, 0, 0, -u4*x4, -u4*y4, u4},
					{0, 0, 0, x4, y4, 1, -v4*x4, -v4*y4, v4},	
				  };
	//if A[0][0] is 0
	if(A[0][0] == 0)
	{
		for(int i = 1;i < 8;i++)
		{
			if(A[i][0] != 0)
			{
				//swap the row and break
				double temp;
				for(int j = 0;j < 9;j++)
				{
					temp = A[0][j];
					A[0][j] = A[i][j];
					A[i][j] = temp;
				}
				break;
			}
		}
	}

	//begin
	for(int i = 1;i < 8;i++)
	{
		// get the high row
		double max = 0;
		int index;
		for(int j = i-1;j < 8;j++)
		{
			if(abs(A[j][i-1]) > max)
			{
				max = abs(A[j][i-1]);
				index = j;
			}
		}
		for(int j = 0;j < 9;j++)
		{
			double temp = A[i-1][j];
			A[i-1][j] = A[index][j];
			A[index][j] = temp;		
		}

		//work on each row
		for(int j = i;j < 8;j++)
		{
			double x = A[j][i-1] / A[i-1][i-1];
			for(int k = i-1;k < 9;k++)
			{
				A[j][k] = A[j][k] - x*A[i-1][k];
			}
		}

		//if A[i][i] = 0
		if(A[i][i] == 0)
		{
			for(int j = i+1;j < 8;j++)
			{
				if(A[j][i] != 0)
				{
					//swap the row and break
					double temp;
					for(int k = 0;k < 9;k++)
					{
						temp = A[i][k];
						A[i][k] = A[j][k];
						A[j][k] = temp;
					}
					break;
				}
			}
		}
	}

	vector<double> ans(8);
	for(int i = 7;i >= 0;i--)
	{
		double b = A[i][8];
		for(int j = 7;j >= i+1;j--)
			b = b - A[i][j] * ans[j];
		ans[i] = b/A[i][i];
	}
	ans.push_back(1);
	return ans;
}

pair<double, double> point_transform(vector<double> para, int u, int v)
{
	double x = para[0]*u + para[1]*v + para[2]*1;
	double y = para[3]*u + para[4]*v + para[5]*1;
	double q = para[6]*u + para[7]*v + para[8]*1;
	x /= q;
	y /= q;
	return pair<double, double>(x, y);
}

vector<double> inverse(vector<double> matrix) 
{
	vector<double> ans;
	double a11 = matrix[0];
	double a12 = matrix[1];
	double a13 = matrix[2];
	double a21 = matrix[3];
	double a22 = matrix[4];
	double a23 = matrix[5];
	double a31 = matrix[6];
	double a32 = matrix[7];
	double a33 = matrix[8];

	double det = a11*a22*a33+a21*a32*a13+a31*a12*a23-a11*a32*a23-a31*a22*a13 - a21*a12*a33;
	if (det != 0) 
	{
		det = 1.0/det;
	}
	ans.push_back(det * (a22*a33-a23*a32) );
	ans.push_back(det * (a13*a32-a12*a33) );
	ans.push_back(det * (a12*a23-a13*a22) );
	ans.push_back(det * (a23*a31-a21*a33) );
	ans.push_back(det * (a11*a33-a13*a31) );
	ans.push_back(det * (a13*a21-a11*a23) );
	ans.push_back(det * (a21*a32-a22*a31) ); 
	ans.push_back(det * (a12*a31-a11*a32) );
	ans.push_back(det * (a11*a22-a12*a21) );
	return ans;
}

vector<match_pair> sift_match(const vector<SiftDescriptor> &desc1, const vector<SiftDescriptor> &desc2)
{
	vector<match_pair> matches;
	for(int i = 0;i < desc1.size();i++)
	{
		match_pair best1;
		match_pair best2;
		double low_dis1 = 10000;
		double low_dis2 = 10000;
		for(int j = 0;j < desc2.size();j++)
		{
			double dis = get_distance(desc1[i].descriptor, desc2[j].descriptor);
			if(dis < low_dis1)
			{
				low_dis1 = dis;
				best1.i = i;
				best1.j = j;
				best1.dis = dis;
			}
			else if(dis < low_dis2)
			{
				low_dis2 = dis;
				best2.i = i;
				best2.j = j;
				best2.dis = dis;
			}
		}
		if(best1.dis < MATCH_THRESH * best2.dis)
		{
			matches.push_back(best1);
		}
	}
	cout << desc1.size() << " " << desc2.size() << endl;
	cout << "Find " << matches.size() << " matches!" << endl;
	return matches;
}

vector<double> ransac(const vector<match_pair> &matches, const vector<SiftDescriptor> &desc1, const vector<SiftDescriptor> &desc2)
{
	vector<double> trans;
	int max_inliner = 0;
	for(int iter = 0;iter < 200;iter++)
	{
		set<pair<double, double> > s;
		vector<int> ran_index;
		while(s.size() < 4)
		{
			int choice = rand() % matches.size();
			double ii = desc1[matches[choice].i].col;
			double jj = desc1[matches[choice].i].row;
			if(s.find(pair<double,double>(ii,jj)) == s.end())
			{
				s.insert(pair<double,double>(ii,jj));
				ran_index.push_back(choice);
			}
		}
		//xy -> uv
		vector<pair<double, double> > uv;
		vector<pair<double, double> > xy;
		for(int i = 0;i < 4;i++)
		{
			int choice = ran_index[i];
			double x = desc1[matches[choice].i].col;
			double y = desc1[matches[choice].i].row;
			double u = desc2[matches[choice].j].col;
			double v = desc2[matches[choice].j].row;
			xy.push_back(pair<double, double>(x,y));
			uv.push_back(pair<double, double>(u,v));
		}
		vector<double> H = get_transform_matrix(xy, uv);

		int inliner = 0;
		for(int k = 0;k < matches.size();k++)
		{
			int x1 = desc1[matches[k].i].col;
			int y1 = desc1[matches[k].i].row;
			int x2 = desc2[matches[k].j].col;
			int y2 = desc2[matches[k].j].row;
			pair<double, double> pt = point_transform(H,x2,y2);
			int x3 = pt.first;
			int y3 = pt.second;
			if( ( (x3-x1)*(x3-x1)+(y3-y1)*(y3-y1) ) <= 30 )
			{
				inliner++;
			}
		}
		if(inliner > max_inliner)
		{
			max_inliner = inliner;
			trans = H;
		}
	}
	cout << "max_inline " << max_inliner << endl;

	cout << trans[0] << " " << trans[1] << " " << trans[2] << endl;
	cout << trans[3] << " " << trans[4] << " " << trans[5] << endl;
	cout << trans[6] << " " << trans[7] << " " << trans[8] << endl;
	return trans;
}

CImg<double> stitch(vector<double> trans, const CImg<double> &img1, const CImg<double> &img2)
{
	vector<double> transs = inverse(trans);
	int width = img1.width();
	int height = img1.height();
	pair<double, double> pt1 = point_transform(trans, 0, 0);
	pair<double, double> pt2 = point_transform(trans, width-1, 0);
	pair<double, double> pt3 = point_transform(trans, 0, height-1);
	pair<double, double> pt4 = point_transform(trans, width-1, height-1);

	int min_x = min(pt1.first, pt3.first);
	int min_y = min(pt1.second, pt2.second);
	int max_x = max(pt2.first, pt4.first);
	int max_y = max(pt3.second, pt4.second);
	cout << min_x << " " << max_x << " " << min_y << " " << max_y << endl;


	int new_width = width - min_x;
	int new_height = max(max_y,height)-min(0,min_y);

	int x_offset = -min_x;
	int y_offset;
	if(min_y < 0)
		y_offset = -min_y;
	else
		y_offset = 0;

	cout << new_width << " " << new_height << endl;
	CImg<double> out1(new_width, new_height, 1, 3, 0);
	CImg<double> out2(new_width, new_height, 1, 3, 0);

	cimg_forXY(out2, x, y)
	{
		pair<double, double> pt = point_transform(transs, x+min_x, y+min_y);
		int xx = pt.first;
		int yy = pt.second;
		if(0 <= xx && xx < img2.width() && 0 <= yy && yy < img2.height())
		{
			out2(x,y,0) = interpolation(img2,xx,yy,0);
			out2(x,y,1) = interpolation(img2,xx,yy,1);
			out2(x,y,2) = interpolation(img2,xx,yy,2);
		}
	}

	cimg_forXY(img1, x, y)
	{
		double r = img1(x,y,0);
		double g = img1(x,y,1);
		double b = img1(x,y,2);
		if(r+g+b > 10)
		{
			out1(x+x_offset,y+y_offset,0) = r;
			out1(x+x_offset,y+y_offset,1) = g;
			out1(x+x_offset,y+y_offset,2) = b;
		}
	}

	CImg<double> out = blendTwoImages(out1, out2);

	width = out.width();
	height = out.height();

	int up = 0;
	int down = height-1;
	while(1)
	{
		bool flag = true;
		for(int x = 0;x < width;x++)
		{
			if(!(out(x,up,0)+out(x,up,1)+out(x,up,2) < 30))
			{
				flag = false;
				break;
			}
		}
		if(!flag)
			break;
		up++;
	}
	while(1)
	{
		bool flag = true;
		for(int x = 0;x < width;x++)
		{
			if(!(out(x,down,0)+out(x,down,1)+out(x,down,2) < 30))
			{
				flag = false;
				break;
			}
		}
		if(!flag)
			break;
		down--;
	}
	out.crop(0,up,0,0,width-1,down,0,2);
	
	return out;
}

CImg<double> img_stitch(const CImg<double> &img1, const CImg<double> &img2)
{
	CImg<double> gray1 = get_gray_scale(img1);
	CImg<double> gray2 = get_gray_scale(img2);
	vector<SiftDescriptor> desc1 = Sift::compute_sift(gray1.crop(0,0,0,0,383,gray1.height()-1,0,0));
	vector<SiftDescriptor> desc2 = Sift::compute_sift(gray2);
	vector<match_pair> matches = sift_match(desc1, desc2);
	vector<double> trans = ransac(matches, desc1, desc2);
	CImg<double> result = stitch(trans, img1, img2);
	cout << "Done!" << endl << endl;
	return result;
}

bool isEmpty(const CImg<unsigned char> &img, int x, int y) 
{
	return (img(x, y, 0) == 0 && img(x, y, 1) == 0 && img(x, y, 2) == 0);
}

CImg<double> blendTwoImages(const CImg<double> &a, const CImg<double> &b) 
{
	double sum_a_x = 0;
	int a_n = 0;

	double sum_overlap_x = 0;
	int overlap_n = 0;

	for (int x = 0; x < a.width(); x++) 
	{
		if (!isEmpty(a, x, a.height() / 2)) 
		{
			sum_a_x += x;
			a_n++;
		}

		if (!isEmpty(a, x, a.height() / 2) && !isEmpty(b, x, a.height() / 2)) 
		{
			sum_overlap_x += x;
			overlap_n++;
		}
	}

	int n_level = floor(log2(a.height()));

	vector<CImg<double> > a_pyramid(n_level);
	vector<CImg<double> > b_pyramid(n_level);
	vector<CImg<double> > mask(n_level);

	a_pyramid[0] = a;
	b_pyramid[0] = b;
	mask[0] = CImg<double>(a.width(), a.height(), 1, 1, 0);

	if (sum_a_x / a_n < sum_overlap_x / overlap_n) 
	{
		for (int x = 0; x < sum_overlap_x / overlap_n; x++) 
		{
			for (int y = 0; y < a.height(); y++) 
			{
				mask[0](x, y) = 1;
			}
		}
	}
	else 
	{
		for (int x = sum_overlap_x / overlap_n + 1; x < a.width(); x++) 
		{
			for (int y = 0; y < a.height(); y++)
			{
				mask[0](x, y) = 1;
			}
		}
	}

	for (int i = 1; i < n_level; i++) 
	{
		a_pyramid[i] = a_pyramid[i - 1].get_blur(2).get_resize(a_pyramid[i - 1].width() / 2, a_pyramid[i - 1].height() / 2, 1, a_pyramid[i - 1].spectrum(), 3);
		b_pyramid[i] = b_pyramid[i - 1].get_blur(2).get_resize(b_pyramid[i - 1].width() / 2, b_pyramid[i - 1].height() / 2, 1, b_pyramid[i - 1].spectrum(), 3);
		
		mask[i] = mask[i - 1].get_blur(2).get_resize(mask[i - 1].width() / 2, mask[i - 1].height() / 2, 1, mask[i - 1].spectrum(), 3);
	}

	for (int i = 0; i < n_level - 1; i++) 
	{
		a_pyramid[i] = a_pyramid[i] - a_pyramid[i + 1].get_resize(a_pyramid[i].width(), a_pyramid[i].height(), 1, a_pyramid[i].spectrum(), 3);
		b_pyramid[i] = b_pyramid[i] - b_pyramid[i + 1].get_resize(b_pyramid[i].width(), b_pyramid[i].height(), 1, b_pyramid[i].spectrum(), 3);
	}

	vector<CImg<double> > blend_pyramid(n_level);

	for (int i = 0; i < n_level; i++) 
	{
		blend_pyramid[i] = CImg<double>(a_pyramid[i].width(), a_pyramid[i].height(), 1, a_pyramid[i].spectrum(), 0);
		for (int x = 0; x < blend_pyramid[i].width(); x++) 
		{
			for (int y = 0; y < blend_pyramid[i].height(); y++) 
			{
				for (int k = 0; k < blend_pyramid[i].spectrum(); k++) 
				{
					blend_pyramid[i](x, y, k) = a_pyramid[i](x, y, k) * mask[i](x, y) + b_pyramid[i](x, y, k) * (1.0 - mask[i](x, y));
				}
			}
		}
	}

	CImg<double> res = blend_pyramid[n_level - 1];
	for (int i = n_level - 2; i >= 0; i--) {
		res.resize(blend_pyramid[i].width(), blend_pyramid[i].height(), 1, blend_pyramid[i].spectrum(), 3);

		for (int x = 0; x < blend_pyramid[i].width(); x++) 
		{
			for (int y = 0; y < blend_pyramid[i].height(); y++) 
			{
				for (int k = 0; k < blend_pyramid[i].spectrum(); k++) 
				{
					double t = res(x, y, k) + blend_pyramid[i](x, y, k);

					if (t > 255) 
					{
						t = 255;
					}
					else if (t < 0) 
					{
						t = 0;
					}

					res(x, y, k) = t;
				}
			}
		}
	}
	return res;
}
