#ifndef _CANNY_H_
#define _CANNY_H_

#include "CImg.h"
using namespace cimg_library;
using namespace std;

class canny{
private:
	CImg<unsigned char> src;   // the source image
	CImg<unsigned char> gray;  // the gray image
    CImg<unsigned char> edge;  // the edge image
    int width;
    int height;

    //code0 struct of canny
    unsigned char *data; 	/* input image */
    int *idata;          	/* output for edges */
    int *magnitude;      	/* edge magnitude as detected by Gaussians */
    float *xConv;        	/* temporary for convolution in x direction */
    float *yConv;        	/* temporary for convolution in y direction */
    float *xGradient;    	/* gradients in x direction, as detected by Gaussians */
    float *yGradient;    	/* gradients in y direction, as detected by Gaussians */

    //parameter
    float lowThreshold;
    float highThreshold;
    int gaussianKernelWidth;
    float gaussianKernelRadius;
    bool contrastNormalised;

public:
	canny(CImg<unsigned char> img);
    canny(CImg<unsigned char> img, int w, int r, int l, int h);
	void rgb2gray();
	void edgedetect();
    void allocatebuffers();
    void killbuffers();
    void normalizeContrast();
    bool computeGradients();
    void performHysteresis();
    void follow(int x1, int y1, int i1, int threshold);
    void showEdgeDetected();
    CImg<unsigned char> getSrc();
    CImg<unsigned char> getEdge();
};


#endif