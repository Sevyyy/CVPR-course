#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include "can.h"
using namespace std;

int main(int argc, char* argv[])
{
	//define the parameter for the 6 image {kernelwidth, radius, low, high, threshold}
	int p0[5] = {16, 3, 2, 4, 550};
	int p1[5] = {7,  3, 4, 8, 510};
	int p2[5] = {7,  3, 4, 8, 550};
	int p3[5] = {5,  2, 2, 4, 550};
	int p4[5] = {5,  2, 2, 4, 550};
	int p5[5] = {5,  3, 2, 4, 400};
	vector< vector<int> > para(6);
	para[0] = vector<int>(p0, p0+5);
	para[1] = vector<int>(p1, p1+5);
	para[2] = vector<int>(p2, p2+5);
	para[3] = vector<int>(p3, p3+5);
	para[4] = vector<int>(p4, p4+5);
	para[5] = vector<int>(p5, p5+5);

	//output file output.txt for each line function in each image
	ofstream fout;
	fout.open("line.txt");
	for(int img = 0;img < 6;img++)
	{
		//calculate some parameter, using canny to get the edge
		string data_path = "./Dataset/";
		int low = (int)(para[img][2]*MAGNITUDE_SCALE+0.5f);
		int high = (int)(para[img][3]*MAGNITUDE_SCALE+0.5f);
		canny can(data_path+"0" + (char)(img+'1') + ".bmp");
		can.to_gray_scale();
		int kwidth = can.gaussian_filter(para[img][0],para[img][1]);
		can.non_maximal_supression(kwidth);
		can.perform_hysteresis(low,high);

		//begin hough
		CImg<unsigned char> src;
		src.load_bmp((data_path+"0" + (char)(img+'1') + ".bmp").c_str());
		CImg<int> edge = can.get_edge_img();
		int width = edge._width;
		int height = edge._height;
		int center_x = width/2;
		int center_y = height/2;
		int num_angle = 360;
		int max_rho = (int)sqrt(center_x*center_x + center_y*center_y);
		vector< vector<int> > accum(num_angle, vector<int>(2*max_rho+1,0));
		vector<double> sinv(num_angle,0);
		vector<double> cosv(num_angle,0);
		for(int i = 0;i < num_angle; i++)
		{
			sinv[i] = sin(i*1.0*cimg::PI / num_angle);
			cosv[i] = cos(i*1.0*cimg::PI / num_angle);
		}

		//transform to hough space
		cimg_forXY(edge, x, y)
		{
			if(edge(x,y) > 0)
			{
				for(int i = 0;i < num_angle;i++)
				{
					int rho = (int)((x-center_x)*cosv[i] + (y-center_y)*sinv[i]);
					rho += max_rho;
					accum[i][rho]++;
				}
			}
		}

		//performing threshold
		cout << para[img][4] << endl;
		CImg<int> hough_space(num_angle, 2*max_rho+1, 1, 1, 0);
		vector< vector<int> > buf;
		for(int i = 0;i < num_angle;i++)
		{
			for(int j = 0;j < 2*max_rho+1;j++)
			{
				
				if(accum[i][j] > para[img][4])
				{
					vector<int> temp;
					temp.push_back(i);
					temp.push_back(j);
					buf.push_back(temp);
				}
				hough_space(i,j) = 5*accum[i][j];
			}
		}
		hough_space.display();
		
		//find local maximum
		CImg<unsigned char> out(width, height,1 ,1 ,0);
		int window_size = 5;
		vector< vector<int> > node;
		fout << "Image" << img+1 << ":" << endl;
		for(int i = 0;i < buf.size();i++)
		{
			int xita = buf[i][0];
			int rho = buf[i][1] - max_rho;
			bool valid = true;
			//create window to find local maximum
			for(int p = -(window_size/2);p < window_size/2+1;p++)
			{
				for(int q = -5*(window_size/2);q < 5*window_size/2+1;q++)
				{
					if(p!=0 || q!=0)
					{
						int x = xita + p;
						int y = rho + max_rho + q;
						if(y < 0)
							continue;
						if(y < 2*max_rho)
						{
							if(x < 0)
							{
								x += num_angle;
								y = -rho+max_rho+q;
							}
							if(x >= num_angle)
							{
								x -= num_angle;
								y = -rho+max_rho+q;
							}
							if(accum[x][y] <= accum[xita][rho+max_rho])
								continue;
							valid = false;
							break;
						}
					}
				}
			}
			if(!valid)
				continue;

			//draw the valid line and save the function
			if(xita >= num_angle/4 && xita <= 3*num_angle/4)
			{
				for(int j = 0;j < width;j++)
				{
					int temp = (rho - (j-center_x)*cosv[xita]) / sinv[xita] + center_y;
					if(temp >= 0 && temp < height)
					{
						if(out(j,temp) == 255)
						{
							vector<int> nd;
							nd.push_back(j);
							nd.push_back(temp);
							node.push_back(nd);
						}
						out(j, temp) = 255;
						src(j, temp, 0) = 255;
						src(j, temp, 1) = 0;
						src(j, temp, 2) = 0;
					}
				}
			}
			else
			{
				for(int j = 0;j < height;j++)
				{
					int temp = (rho - (j-center_y)*sinv[xita]) / cosv[xita] + center_x;
					if(temp >= 0 && temp < width)	
					{
						if(out(temp, j) == 255)
						{
							vector<int> nd;
							nd.push_back(temp);
							nd.push_back(j);
							node.push_back(nd);
						}
						out(temp, j) = 255;
						src(temp, j, 0) = 255;
						src(temp, j, 1) = 0;
						src(temp, j, 2) = 0;
					}
				}
			}

			cout << "xita : " << xita << " rho : " << rho << " : " << accum[xita][rho+max_rho] <<endl;
			cout << "Line : " << cosv[xita] << "*x + " << sinv[xita] << "*y = " << rho + center_x*cosv[xita] + center_y*sinv[xita] << endl;
			fout << "  Line : " << cosv[xita] << "*x + " << sinv[xita] << "*y = " << rho + center_x*cosv[xita] + center_y*sinv[xita] << endl;
		}

		const unsigned char red[] = {255,0,0};
		for(int i = 0;i < node.size();i++)
		{
			src.draw_circle(node[i][0], node[i][1], 40, red);
		}

		(src,edge,out).display();
		string out_path = "./out/";
		src.save((out_path + (char)(img+'1') + "_out.bmp").c_str());
	}
	fout.close();
	return 0;
}