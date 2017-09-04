#include "CImg.h"
#include "canny.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace std;
using namespace cimg_library;

//define some paremeter
const int threshold = 300;
const int num_angle = 360;
const unsigned char pink[] = {241,158,194};
//the num_angle values of the sin and cos, avoid computing repeatedly
vector<double> sinv(num_angle,0);
vector<double> cosv(num_angle,0);

// for transform the index to filename
string int_to_string(int n);

//compare function of the type pair<pair<int,int>, int>, used in sort()
bool comp(const pair<pair<int,int>, int> &a, const pair<pair<int,int>, int> &b);

//get the 3*3matrix by Gaussian elimination on the 8 points(8 equation)
vector<double> get_transform_matrix(vector<pair<double, double> > uv, vector<pair<double,double> > xy);

//multiply the (u,v) and A to get (x,y), another word is map
pair<double, double> mul(vector<double> para, int u, int v);

//transform to gray scale and get the edge image with blur
CImg<float> get_edge(CImg<unsigned char> src);

//transform to hough space
vector<vector<int> > to_hough(CImg<float> edge);

//get local maximum
vector<pair<pair<int,int>, int> > get_local_maximum(vector<vector<int> > accum);

//draw lines and points and get the nodes
vector<pair<double, double> > draw_lines_and_points(vector<pair<pair<int,int>, int> > lines, CImg<unsigned char> src, CImg<unsigned char> &src_out);

//get the a4 paper and display the result and return the a4 paper
CImg<unsigned char> get_a4_image(CImg<unsigned char> src, vector<pair<double, double> > node, CImg<unsigned char> src_out);

int main()
{
	//the num_angle values of the sin and cos, avoid computing repeatedly
	for(int i = 0;i < num_angle; i++)
	{
		sinv[i] = sin(i*2.0*cimg::PI / num_angle);
		cosv[i] = cos(i*2.0*cimg::PI / num_angle);
	}

	//begin the loop of the 16 image
	for(int image = 1;image <= 16;image++)
	{
		//load file and do edge detect
		string filepath = "./Dataset/" + int_to_string(image) + ".bmp";
		CImg<unsigned char> src(filepath.c_str());
		CImg<unsigned char> src_out = src;                  //for image output
		CImg<float> edge = get_edge(src);

		//transform to hough space
		vector< vector<int> > accum = to_hough(edge);

		//find local maximum
		vector<pair<pair<int,int>, int> > lines = get_local_maximum(accum);

		//draw lines and points and get the nodes
		vector<pair<double, double> > node = draw_lines_and_points(lines, src, src_out);

		//get the a4 paper and display the result and save the a4 paper
		CImg<unsigned char> paper =  get_a4_image(src, node, src_out);

		paper.save(("./out/" + int_to_string(image) + "_out.bmp").c_str());
		(src_out,paper).display();
	}
	return 0;
}

// for transform the index to filename
string int_to_string(int n)
{
	string ans = "";
	while(n)
	{
		ans = (char)((n%10)+'0') + ans;
		n /= 10;
	}
	return ans;
}

//compare function of the type pair<pair<int,int>, int>, used in sort()
bool comp(const pair<pair<int,int>, int> &a, const pair<pair<int,int>, int> &b)
{
	return a.second > b.second;
}

//get the 3*3matrix by Gaussian elimination on the 8 points(8 equation)
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
	return ans;
}

//multiply the (u,v) and A to get (x,y), another word is map
pair<double, double> mul(vector<double> para, int u, int v)
{
	double x = para[0]*u + para[1]*v + para[2]*1;
	double y = para[3]*u + para[4]*v + para[5]*1;
	double q = para[6]*u + para[7]*v + 1;
	x /= q;
	y /= q;
	return pair<double, double>(x, y);
}

//transform to gray scale and get the edge image with blur
CImg<float> get_edge(CImg<unsigned char> src)
{
	CImg<float> edge(src._width, src._height, 1, 1, 0);
	CImg<unsigned char> gray(src._width, src._height, 1, 1, 0);
	cimg_forXY(src, x, y)
	{
		gray(x,y) = 0.2116 * src(x, y, 0) + 0.7152 * src(x, y, 1) + 0.0722 * src(x, y, 2);
	}
	gray.blur(3);
	CImg_3x3(I, float);
	cimg_for3x3(gray, x, y, 0, 0, I, float)
	{
		double ix = Inc - Ipc;
		double iy = Icp - Icn;
		double grad = sqrt(ix*ix + iy*iy);
		if(grad > 15)
			edge(x,y) = 255;
	}
	return edge;
}

//transform to hough space
vector<vector<int> > to_hough(CImg<float> edge)
{
	int width = edge._width;
	int height = edge._height;
	int max_rho = (int)sqrt(width*width + height*height);
	vector< vector<int> > accum(num_angle, vector<int>(max_rho,0));
	
	//transform to hough space
	cimg_forXY(edge, x, y)
	{
		if(edge(x,y) > 0)
		{
			for(int i = 0;i < num_angle;i++)
			{
				int rho = (int)(x*cosv[i] + y*sinv[i]);
				if(rho > 0)
					accum[i][rho]++;
			}
		}
	}
	return accum;
}

//get local maximum
vector<pair<pair<int,int>, int> > get_local_maximum(vector<vector<int> > accum)
{
	int max_rho = accum[0].size();
	//performing threshold
	vector<pair<pair<int,int>, int> > buf; //((x,y),weight)
	for(int i = 0;i < num_angle;i++)
	{
		for(int j = 0;j < max_rho;j++)
		{
			if(accum[i][j] > threshold)
			{
				pair<int,int> temp(i,j);
				buf.push_back(pair<pair<int,int>, int>(temp, accum[i][j]));
			}
		}
	}

	//get local maximum
	int xita_window = 60;
	int rho_windos = 400;
	vector<pair<pair<int,int>, int> > lines;
	for(int i = 0;i < buf.size();i++)
	{
		int xita = buf[i].first.first;
		int rho = buf[i].first.second;
		bool valid = true;
		//create window to find local maximum
		for(int p = -xita_window;p < xita_window+1;p++)
		{
			for(int q = -rho_windos;q < rho_windos+1;q++)
			{
				if(p!=0 || q!=0)
				{
					int x = xita + p;
					int y = rho + q;
					if(y < 0)
						continue;
					if(y < max_rho)
					{
						if(x < 0)
						{
							x += num_angle;
						}
						if(x >= num_angle)
						{
							x -= num_angle;
						}
						if(accum[x][y] < accum[xita][rho])
							continue;
						if(accum[x][y] == accum[xita][rho])
						{
							accum[xita][rho]++;
							continue;
						}
						valid = false;
						break;
					}
				}
			}
			if(!valid)
				break;
		}
		if(!valid)
			continue;
		//recode the line detected
		lines.push_back(buf[i]);
	}
	
	//sort it
	sort(lines.begin(), lines.end(), comp);
	return lines;
}

//draw lines and points and get the nodes
vector<pair<double, double> > draw_lines_and_points(vector<pair<pair<int,int>, int> > lines, CImg<unsigned char> src, CImg<unsigned char> &src_out)
{
	int width = src._width;
	int height = src._height;
	int max_rho = (int)sqrt(width*width + height*height);
	//for lines
	for(int i = 0;i < 4;i++)
	{
		int xita = lines[i].first.first;
		int rho = lines[i].first.second;
		double a = cosv[xita];
		double b = sinv[xita];
		double c = rho;
		vector<pair<int, int> > point;
		if(a == 0)
		{
			point.push_back(pair<int, int>(0,c/b));
			point.push_back(pair<int, int>(width-1,c/b));
		}
		else if(b == 0)
		{
			point.push_back(pair<int, int>(c/a,0));
			point.push_back(pair<int, int>(c/a,height-1));
		}
		else
		{
			if(0 < c/a && c/a < width-1)
			{
				point.push_back(pair<int, int>(c/a,0));
			}
			if(0 < (c-b*(height-1))/a && (c-b*(height-1))/a < width-1)
			{
				point.push_back(pair<int, int>((c-b*(height-1))/a,height-1));
			}
			if(0 < c/b && c/b < height-1)
			{
				point.push_back(pair<int, int>(0,c/b));
			}
			if(0 < (c-a*(width-1))/b && (c-a*(width-1))/b < height-1)
			{
				point.push_back(pair<int, int>(width-1,(c-a*(width-1))/b));
			}
		}
		src_out.draw_line(point[0].first, point[0].second, point[1].first, point[1].second, pink);
	}

	//compute node position
	vector<pair<double, double> > node;
	for(int i = 0;i < 4;i++)
	{
		for(int j = i+1;j < 4;j++)
		{
			int xitai = lines[i].first.first;
			int rhoi = lines[i].first.second;
			int xitaj = lines[j].first.first;
			int rhoj = lines[j].first.second;
			double a1 = cosv[xitai];
			double b1 = sinv[xitai];
			double c1 = rhoi;
			double a2 = cosv[xitaj];
			double b2 = sinv[xitaj];
			double c2 = rhoj;
			int y = (int)((a1*c2-a2*c1)/(a1*b2-a2*b1));
			int x = (int)((b1*c2-b2*c1)/(b1*a2-b2*a1));
			if(x >= 0 && x <= width && y >= 0 && y <= height)
			{
				node.push_back(pair<double, double>(x,y));
			}
		}
	}
	for(int i = 0;i < node.size();i++)
	{
		src_out.draw_circle((int)(node[i].first), (int)(node[i].second), 40, pink);
	}

	return node;
}

//get the a4 paper and display the result and return the a4 paper
CImg<unsigned char> get_a4_image(CImg<unsigned char> src, vector<pair<double, double> > node, CImg<unsigned char> src_out)
{
	int width = src._width;
	int height = src._height;
	int max_rho = (int)sqrt(width*width + height*height);

	vector<pair<double,double> > uv;
	CImg<unsigned char> paper;
	//vertical
	if(width < height)
	{
		uv.push_back(pair<double, double>(0,0));
		uv.push_back(pair<double, double>(1999,0));
		uv.push_back(pair<double, double>(0,2827));
		uv.push_back(pair<double, double>(1999,2827));
		paper = CImg<unsigned char>(2000,2828,1,3,0);
	}
	//horizontal
	else
	{
		uv.push_back(pair<double, double>(0,0));
		uv.push_back(pair<double, double>(2827,0));
		uv.push_back(pair<double, double>(0,1999));
		uv.push_back(pair<double, double>(2827,1999));
		paper = CImg<unsigned char>(2828,2000,1,3,0);
	}

	//match the 4 point
	vector< pair<double, double> > xy(4);
	for(int i = 0;i < 4;i++)
	{
		int flag1 = 0;
		int flag2 = 0;
		for(int j = 0;j < 4;j++)
		{
			if(i != j)
			{
				if(node[i].first >= node[j].first)
					flag1++;
				if(node[i].second >= node[j].second)
					flag2++;
			}
		}
		if(flag1 >= 2)
		{
			if(flag2 < 2)
				xy[1] = node[i];
			else
				xy[3] = node[i];
		}
		else
		{
			if(flag2 < 2)
				xy[0] = node[i];
			else
				xy[2] = node[i];
		}
	}
	vector<double> parameter = get_transform_matrix(xy, uv);
	cimg_forXY(paper, x, y)
	{
		pair<double, double> point = mul(parameter, x, y);
		double x0 = point.first;
		double y0 = point.second;
		int px = (int)x0;
		int py = (int)y0;
		double a = x0-px;
		double b = y0-py;
		paper(x,y,0) = (unsigned char)((1-a)*(1-b)*src(px,py,0) + a*(1-b)*src(px+1,py,0) + (1-a)*b*src(px,py+1,0) + a*b*src(px+1,py+1,0));
		paper(x,y,1) = (unsigned char)((1-a)*(1-b)*src(px,py,1) + a*(1-b)*src(px+1,py,1) + (1-a)*b*src(px,py+1,1) + a*b*src(px+1,py+1,1));
		paper(x,y,2) = (unsigned char)((1-a)*(1-b)*src(px,py,2) + a*(1-b)*src(px+1,py,2) + (1-a)*b*src(px,py+1,2) + a*b*src(px+1,py+1,2));
	}
	return paper;
}