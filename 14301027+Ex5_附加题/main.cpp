#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstring>
#include "CImg.h"

using namespace std;
using namespace cimg_library;

//define the triangle class to represent a triangle by three points
class triangle
{
public:
	int x0, x1, x2, y0, y1, y2;
	//constructor
	triangle(int x0, int y0, int x1, int y1, int x2, int y2);
	//judge if a point x,y is in the triangle
	bool contain(int x, int y);
	//print out the three points
	void print();
};

//load the triangle from txt
void load_data(vector<triangle> &ori, vector<triangle> &tar);

//get the affine transfrom parameter
vector<double> get_transform_parameter(triangle ori, triangle tar);

//do transform for a point
pair<int, int> do_transform(vector<double> para, int x, int y);

// for transform the index to filename
string int_to_string(int n);


int main()
{
	//load the triangle data from triangle.txt
	vector<triangle> ori;          //the origin image
	vector<triangle> tar;          //the target image
	load_data(ori, tar);

	//load the two image
	CImg<unsigned char> man("origin.bmp");
	CImg<unsigned char> girl("target.bmp");

	vector<CImg<unsigned char> > frame;
	vector<double> p;
	for(int i = 1;i <= 11;i++)
	{
		p.push_back(i*1.0/12);
	}

	for(int t = 0;t < 11;t++)
	{
		CImg<unsigned char> temp(man._width, man._height, 1, 3, 0);

		//get the medium triangle by p
		vector<triangle> med;
		for(int i = 0;i < ori.size();i++)
		{
			int x0 = (int)((1-p[t])*ori[i].x0 + p[t]*tar[i].x0);
			int y0 = (int)((1-p[t])*ori[i].y0 + p[t]*tar[i].y0);
			int x1 = (int)((1-p[t])*ori[i].x1 + p[t]*tar[i].x1);
			int y1 = (int)((1-p[t])*ori[i].y1 + p[t]*tar[i].y1);
			int x2 = (int)((1-p[t])*ori[i].x2 + p[t]*tar[i].x2);
			int y2 = (int)((1-p[t])*ori[i].y2 + p[t]*tar[i].y2);
			med.push_back(triangle(x0, y0, x1, y1, x2, y2));
		}

		//get the parameter for transformation from med to ori and tar
		vector<vector<double> > para_ori;
		vector<vector<double> > para_tar;
		for(int i = 0;i < ori.size();i++)
		{
			para_ori.push_back(get_transform_parameter(med[i], ori[i]));
			para_tar.push_back(get_transform_parameter(med[i], tar[i]));
		}

		//get the medium iamge
		cimg_forXY(temp, x, y)
		{
			for(int i = 0;i < med.size();i++)
			{
				if(med[i].contain(x,y))
				{
					pair<int, int> pos_ori = do_transform(para_ori[i], x, y);
					pair<int, int> pos_tar = do_transform(para_tar[i], x, y);
					temp(x,y,0) = (unsigned char)( (1-p[t])*man(pos_ori.first,pos_ori.second,0) + p[t]*girl(pos_tar.first,pos_tar.second,0) );
					temp(x,y,1) = (unsigned char)( (1-p[t])*man(pos_ori.first,pos_ori.second,1) + p[t]*girl(pos_tar.first,pos_tar.second,1) );
					temp(x,y,2) = (unsigned char)( (1-p[t])*man(pos_ori.first,pos_ori.second,2) + p[t]*girl(pos_tar.first,pos_tar.second,2) );
					break;
				}
			}
		}
		temp.display();
		temp.save(("./out/" + int_to_string(t+1) + ".bmp").c_str());
	}
	return 0;
}

//constructor
triangle::triangle(int x0, int y0, int x1, int y1, int x2, int y2)
{
	this->x0 = x0;
	this->y0 = y0;
	this->x1 = x1;
	this->y1 = y1;
	this->x2 = x2;
	this->y2 = y2;
}

//judge if a point x,y is in the triangle
bool triangle::contain(int x, int y)
{
	int a1 = x1-x0;
	int a2 = y1-y0;
	int b1 = x2-x0;
	int b2 = y2-y0;
	int c1 = x -x0;
	int c2 = y -y0;
	double u = (c1*b2-c2*b1)*1.0/(a1*b2-a2*b1);
	double v = (a1*c2-a2*c1)*1.0/(a1*b2-a2*b1);
	if(u >= 0 && v >= 0 && u+v <= 1)
		return true;
	else
		return false;
}

//print out the three points
void triangle::print()
{
	cout << x0 << ", " << y0 << ", " << x1 << ", " << y1 << ", " << x2 << ", " << y2 << endl;
	return;
}

//get the affine transfrom parameter
vector<double> get_transform_parameter(triangle ori, triangle tar)
{
	vector<double> ret(6);
	ret[0] = 1.0 * (tar.x0*ori.y1 - tar.x1*ori.y0 - tar.x0*ori.y2 + tar.x2*ori.y0 + tar.x1*ori.y2 - tar.x2*ori.y1)/(ori.x0*ori.y1 - ori.x1*ori.y0 - ori.x0*ori.y2 + ori.x2*ori.y0 + ori.x1*ori.y2 - ori.x2*ori.y1);
	ret[1] = -1.0 * (tar.x0*ori.x1 - tar.x1*ori.x0 - tar.x0*ori.x2 + tar.x2*ori.x0 + tar.x1*ori.x2 - tar.x2*ori.x1)/(ori.x0*ori.y1 - ori.x1*ori.y0 - ori.x0*ori.y2 + ori.x2*ori.y0 + ori.x1*ori.y2 - ori.x2*ori.y1);
	ret[2] = 1.0 * (tar.x0*ori.x1*ori.y2 - tar.x0*ori.x2*ori.y1 - tar.x1*ori.x0*ori.y2 + tar.x1*ori.x2*ori.y0 + tar.x2*ori.x0*ori.y1 - tar.x2*ori.x1*ori.y0)/(ori.x0*ori.y1 - ori.x1*ori.y0 - ori.x0*ori.y2 + ori.x2*ori.y0 + ori.x1*ori.y2 - ori.x2*ori.y1);
	ret[3] = 1.0 * (tar.y0*ori.y1 - tar.y1*ori.y0 - tar.y0*ori.y2 + tar.y2*ori.y0 + tar.y1*ori.y2 - tar.y2*ori.y1)/(ori.x0*ori.y1 - ori.x1*ori.y0 - ori.x0*ori.y2 + ori.x2*ori.y0 + ori.x1*ori.y2 - ori.x2*ori.y1);
	ret[4] = -1.0 * (tar.y0*ori.x1 - tar.y1*ori.x0 - tar.y0*ori.x2 + tar.y2*ori.x0 + tar.y1*ori.x2 - tar.y2*ori.x1)/(ori.x0*ori.y1 - ori.x1*ori.y0 - ori.x0*ori.y2 + ori.x2*ori.y0 + ori.x1*ori.y2 - ori.x2*ori.y1);
	ret[5] = 1.0 * (tar.y0*ori.x1*ori.y2 - tar.y0*ori.x2*ori.y1 - tar.y1*ori.x0*ori.y2 + tar.y1*ori.x2*ori.y0 + tar.y2*ori.x0*ori.y1 - tar.y2*ori.x1*ori.y0)/(ori.x0*ori.y1 - ori.x1*ori.y0 - ori.x0*ori.y2 + ori.x2*ori.y0 + ori.x1*ori.y2 - ori.x2*ori.y1);
	return ret;
}

//load the triangle from txt
void load_data(vector<triangle> &ori, vector<triangle> &tar)
{
	vector<pair<int, int> > p1;
	vector<pair<int, int> > p2;
	ifstream fin;
	fin.open("point_ori.txt");
	int x,y,z;
	while(fin >> x)
	{
		fin >> y;
		p1.push_back(pair<int, int>(x,y));
	}
	fin.close();
	fin.open("point_tar.txt");
	while(fin >> x)
	{
		fin >> y;
		p2.push_back(pair<int, int>(x,y));
	}
	fin.close();
	fin.open("triangle.txt");
	while(fin >> x)
	{
		fin >> y >> z;
		ori.push_back(triangle(p1[x-1].first, p1[x-1].second, p1[y-1].first, p1[y-1].second, p1[z-1].first, p1[z-1].second));
		tar.push_back(triangle(p2[x-1].first, p2[x-1].second, p2[y-1].first, p2[y-1].second, p2[z-1].first, p2[z-1].second));
	}
	fin.close();
}

//do transform for a point
pair<int, int> do_transform(vector<double> para, int x, int y)
{
	int u = (int)(para[0]*x+para[1]*y+para[2]);
	int v = (int)(para[3]*x+para[4]*y+para[5]);
	return pair<int,int>(u,v);
}

// for transform the index to filename
string int_to_string(int n)
{
	if(n == 0)
		return "0";
	string ans = "";
	while(n)
	{
		ans = (char)((n%10)+'0') + ans;
		n /= 10;
	}
	return ans;
}