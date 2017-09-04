#include <cstdlib>
#include <iostream>
#include "canny.h"
#include "CImg.h"

using namespace cimg_library;
using namespace std;

#define ffabs(x) ( (x) >= 0 ? (x) : -(x) ) 
#define GAUSSIAN_CUT_OFF 0.005f
#define MAGNITUDE_SCALE 100.0f
#define MAGNITUDE_LIMIT 1000.0f
#define MAGNITUDE_MAX ((int) (MAGNITUDE_SCALE * MAGNITUDE_LIMIT))

float hypotenuse(float x, float y)
{
	return (float) sqrt(x*x + y*y);
}

float gaussian(float x, float sigma)
{
	return (float) exp(-(x*x) / (2.0f*sigma*sigma));
}

canny::canny(CImg<unsigned char> img) 
{
	this->src = img;
	rgb2gray();
	this->width = this->src.width();
	this->height = this->src.height();

	this->lowThreshold = 10.0f;
    this->highThreshold = 20.0f;
    this->gaussianKernelWidth = 16;
    this->gaussianKernelRadius = 2.0f;
    this->contrastNormalised = false;
}

canny::canny(CImg<unsigned char> img, int w, int r, int l, int h)
{
	this->src = img;
	rgb2gray();
	this->width = this->src.width();
	this->height = this->src.height();

	this->lowThreshold = l;
    this->highThreshold = h;
    this->gaussianKernelWidth = w;
    this->gaussianKernelRadius = r;
    this->contrastNormalised = false;
}


void canny::rgb2gray() 
{
	CImg<unsigned char> imgTemp((this->src).width(), (this->src).height(), 1, 1, 0);
	cimg_forXY(this->src, x, y) 
	{
		imgTemp(x, y, 0) = 0.2116 * this->src(x, y, 0) + 0.7152 * this->src(x, y, 1) + 0.0722 * this->src(x, y, 2);
	};
	this->gray = imgTemp;
}

void canny::edgedetect() 
{
	int err;
	int i;

	allocatebuffers();

	if (contrastNormalised) 
		normalizeContrast();

	if (!computeGradients()) {
		killbuffers();
	}


	int low = (int) (this->lowThreshold * MAGNITUDE_SCALE + 0.5f);
	int high = (int) ( this->highThreshold * MAGNITUDE_SCALE + 0.5f);
	performHysteresis();
	for (int i = 0; i < this->width * this->height; i ++)
		this->idata[i] = this->idata[i] > 0 ? 255 : 0;
	this->edge = CImg<unsigned char>(this->idata, (this->gray).width(), (this->gray).height(), 1, 1, false);
	killbuffers();
	//this->edge.display();
}

void canny::allocatebuffers()
{
	this->data = new unsigned char[width * height];
	this->idata = new int[width * height];
	this->magnitude = new int[width * height];
	this->xConv = new float[width * height];
	this->yConv = new float[width * height];
	this->xGradient = new float[width * height];
	this->yGradient = new float[width * height];
	if(!this->data || !this->idata || !this->magnitude || !this->xConv || !this->yConv || 
		!this->xGradient || !this->yGradient)
		killbuffers();

	memcpy(this->data, (this->gray).data(), width * height);
}

void canny::killbuffers()
{
	if (data)
		delete [] data;
	if (idata)
		delete [] idata;
	if (magnitude)
		delete [] magnitude;
	if (xConv)
		delete [] xConv;
	if (yConv)
		delete [] yConv;
	if (xGradient)
		delete [] xGradient;
	if (yGradient)
		delete [] yGradient;
}

bool canny::computeGradients() 
{	
	float *kernel;
	float *diffKernel;
	int kwidth;

	int initX;
	int maxX;
	int initY;
	int maxY;

	int x, y;
	int i;
	int flag;

	kernel = new float[gaussianKernelWidth];
	diffKernel = new float[gaussianKernelWidth];
	if(!kernel || !diffKernel) {
		goto error_exit;
	}


	/* initialise the Gaussian kernel */
	for (kwidth = 0; kwidth < gaussianKernelWidth; kwidth++) 
	{
		float g1 = gaussian((float)kwidth, gaussianKernelRadius);
		if (g1 <= GAUSSIAN_CUT_OFF && kwidth >= 2) 
			break;
		float g2 = gaussian(kwidth - 0.5f, gaussianKernelRadius);
		float g3 = gaussian(kwidth + 0.5f, gaussianKernelRadius);
		kernel[kwidth] = (g1 + g2 + g3) / 3.0f / (2.0f * (float) 3.14 * gaussianKernelRadius * gaussianKernelRadius);
		diffKernel[kwidth] = g3 - g2;
	}

	initX = kwidth - 1;
	maxX = width - (kwidth - 1);
	initY = width * (kwidth - 1);
	maxY = this->width * (this->height - (kwidth - 1));
	
	/* perform convolution in x and y directions */
	for(x = initX; x < maxX; x++) 
	{
		for(y = initY; y < maxY; y += width) 
		{
			int index = x + y;
			float sumX = this->data[index] * kernel[0];
			float sumY = sumX;
			int xOffset = 1;
			int yOffset = this->width;
			while(xOffset < kwidth) 
			{
				sumY += kernel[xOffset] * (this->data[index - yOffset] + this->data[index + yOffset]);
				sumX += kernel[xOffset] * (this->data[index - xOffset] + this->data[index + xOffset]);
				yOffset += this->width;
				xOffset++;
			}
			
			this->yConv[index] = sumY;
			this->xConv[index] = sumX;
		}

	}

	for (x = initX; x < maxX; x++) 
	{
		for (y = initY; y < maxY; y += this->width) 
		{
			float sum = 0.0f;
			int index = x + y;
			for (i = 1; i < kwidth; i++)
				sum += diffKernel[i] * (this->yConv[index - i] - this->yConv[index + i]);

			this->xGradient[index] = sum;
		}

	}

	for(x = kwidth; x < width - kwidth; x++) 
	{
		for (y = initY; y < maxY; y += width) 
		{
			float sum = 0.0f;
			int index = x + y;
			int yOffset = width;
			for (i = 1; i < kwidth; i++) 
			{
				sum += diffKernel[i] * (this->xConv[index - yOffset] - this->xConv[index + yOffset]);
				yOffset += this->width;
			}

			this->yGradient[index] = sum;
		}
	}

	initX = kwidth;
	maxX = width - kwidth;
	initY = width * kwidth;
	maxY = this->width * (this->height - kwidth);
	for(x = initX; x < maxX; x++) 
	{
		for(y = initY; y < maxY; y += width) 
		{
			int index = x + y;
			int indexN = index - width;
			int indexS = index + width;
			int indexW = index - 1;
			int indexE = index + 1;
			int indexNW = indexN - 1;
			int indexNE = indexN + 1;
			int indexSW = indexS - 1;
			int indexSE = indexS + 1;
			
			float xGrad = this->xGradient[index];
			float yGrad = this->yGradient[index];
			float gradMag = hypotenuse(xGrad, yGrad);

			/* perform non-maximal supression */
			float nMag = hypotenuse(this->xGradient[indexN], this->yGradient[indexN]);
			float sMag = hypotenuse(this->xGradient[indexS], this->yGradient[indexS]);
			float wMag = hypotenuse(this->xGradient[indexW], this->yGradient[indexW]);
			float eMag = hypotenuse(this->xGradient[indexE], this->yGradient[indexE]);
			float neMag = hypotenuse(this->xGradient[indexNE], this->yGradient[indexNE]);
			float seMag = hypotenuse(this->xGradient[indexSE], this->yGradient[indexSE]);
			float swMag = hypotenuse(this->xGradient[indexSW], this->yGradient[indexSW]);
			float nwMag = hypotenuse(this->xGradient[indexNW], this->yGradient[indexNW]);
			float tmp;
			/*
			 * An explanation of what's happening here, for those who want
			 * to understand the source: This performs the "non-maximal
			 * supression" phase of the Canny edge detection in which we
			 * need to compare the gradient magnitude to that in the
			 * direction of the gradient; only if the value is a local
			 * maximum do we consider the point as an edge candidate.
			 * 
			 * We need to break the comparison into a number of different
			 * cases depending on the gradient direction so that the
			 * appropriate values can be used. To avoid computing the
			 * gradient direction, we use two simple comparisons: first we
			 * check that the partial derivatives have the same sign (1)
			 * and then we check which is larger (2). As a consequence, we
			 * have reduced the problem to one of four identical cases that
			 * each test the central gradient magnitude against the values at
			 * two points with 'identical support'; what this means is that
			 * the geometry required to accurately interpolate the magnitude
			 * of gradient function at those points has an identical
			 * geometry (upto right-angled-rotation/reflection).
			 * 
			 * When comparing the central gradient to the two interpolated
			 * values, we avoid performing any divisions by multiplying both
			 * sides of each inequality by the greater of the two partial
			 * derivatives. The common comparand is stored in a temporary
			 * variable (3) and reused in the mirror case (4).
			 * 
			 */
			flag = ( (xGrad * yGrad <= 0.0f) /*(1)*/
				? ffabs(xGrad) >= ffabs(yGrad) /*(2)*/
					? (tmp = ffabs(xGrad * gradMag)) >= ffabs(yGrad * neMag - (xGrad + yGrad) * eMag) /*(3)*/
						&& tmp > fabs(yGrad * swMag - (xGrad + yGrad) * wMag) /*(4)*/
					: (tmp = ffabs(yGrad * gradMag)) >= ffabs(xGrad * neMag - (yGrad + xGrad) * nMag) /*(3)*/
						&& tmp > ffabs(xGrad * swMag - (yGrad + xGrad) * sMag) /*(4)*/
				: ffabs(xGrad) >= ffabs(yGrad) /*(2)*/
					? (tmp = ffabs(xGrad * gradMag)) >= ffabs(yGrad * seMag + (xGrad - yGrad) * eMag) /*(3)*/
						&& tmp > ffabs(yGrad * nwMag + (xGrad - yGrad) * wMag) /*(4)*/
					: (tmp = ffabs(yGrad * gradMag)) >= ffabs(xGrad * seMag + (yGrad - xGrad) * sMag) /*(3)*/
						&& tmp > ffabs(xGrad * nwMag + (yGrad - xGrad) * nMag) /*(4)*/
				);
            if(flag)
			{
				this->magnitude[index] = (gradMag >= MAGNITUDE_LIMIT) ? int(MAGNITUDE_SCALE * MAGNITUDE_LIMIT) : (int) (MAGNITUDE_SCALE * gradMag);
				/*NOTE: The orientation of the edge is not employed by this
				 implementation. It is a simple matter to compute it at
				this point as: Math.atan2(yGrad, xGrad); */
			} 
			else 
			{
				this->magnitude[index] = 0;
			}
		}
	}
	delete[] kernel;
	delete[] diffKernel;
	return true;

error_exit:
		delete[] kernel;
		delete[] diffKernel;
		return false;
}

void canny::performHysteresis() 
{
	int offset = 0;
	int x, y;

	memset(this->idata, 0, this->width * this->height * sizeof(int));
		
	for(y = 0; y < this->height; y++)
	{
		for(x = 0; x < this->width; x++)
		{
		  	if(this->idata[offset] == 0 && this->magnitude[offset] >= this->highThreshold * MAGNITUDE_SCALE + 0.5f) 
		    	follow(x, y, offset, this->lowThreshold * MAGNITUDE_SCALE + 0.5f);
		  	offset++;
	}
	}
}

void canny::follow(int x1, int y1, int i1, int threshold) 
{
	int x, y;
	int x0 = x1 == 0 ? x1 : x1 - 1;
	int x2 = x1 == this->width - 1 ? x1 : x1 + 1;
	int y0 = y1 == 0 ? y1 : y1 - 1;
	int y2 = y1 == this->height -1 ? y1 : y1 + 1;
		
	this->idata[i1] = this->magnitude[i1];
	for (x = x0; x <= x2; x++) 
	{
		for (y = y0; y <= y2; y++) 
		{
	  		int i2 = x + y * this->width;
	  		if ((y != y1 || x != x1) && this->idata[i2] == 0 && this->magnitude[i2] >= threshold) 
	    		follow(x, y, i2, threshold);
		}
	}
}

void canny::normalizeContrast() 
{
	int histogram[256] = {0};
    int remap[256];
	int sum = 0;
    int j = 0;
	int k;
    int target;
    int i;

	for (i = 0; i < width * height; i++) 
			histogram[data[i]]++;
		
	
    for (i = 0; i < 256; i++) 
	{
			sum += histogram[i];
			target = (sum*255)/(width * height);
			for (k = j+1; k <= target; k++) 
				remap[k] = i;
			j = target;
	 }
		
    for (i = 0; i < width * height; i++) 
			data[i] = remap[data[i]];
}

CImg<unsigned char> canny::getSrc() 
{
	return this->src;
}

CImg<unsigned char> canny::getEdge() 
{
	return this->edge;
}