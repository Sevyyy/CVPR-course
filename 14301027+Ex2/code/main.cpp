#include <vector>
#include <algorithm>
#include <string>
#include <cstring>
#include <iostream>
#include "can.h"
using namespace std;

int main(int argc, char* argv[])
{
	string file_path = "../test_Data/";
	string image[4] = { "bigben", "lena", "stpietro", "twows"};
	int i;
	int a, b;
	float c, d;
	cin >> i >> a >> b >> c >> d;
		canny img(file_path + image[i] + ".bmp", image[i]);
		img.to_gray_scale();
		int kwidth = img.gaussian_filter(a,b);
		img.non_maximal_supression(kwidth);
		int low = (int)(c * MAGNITUDE_SCALE + 0.5f);
		int high = (int)(d * MAGNITUDE_SCALE + 0.5f);
		img.perform_hysteresis(low, high);
		cout << "Done!" << i << endl;
	
	return 0;
}