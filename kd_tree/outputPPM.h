#ifndef OUTPUTPPM
#define OUTPUTPPM


#include <iostream>
#include <assert.h>
#include <stdio.h>

//-----------------------create ppm image functon Start--------------------------------------------
struct Pixel {
	unsigned char R, G, B;  // Blue, Green, Red
};

class ColorImage {
	Pixel *pPixel;
	int xRes, yRes;
public:
	ColorImage();  //initialize
	~ColorImage(); //initialize
	void init(int xSize, int ySize);
	void clear(Pixel background);
	Pixel readPixel(int x, int y);
	void writePixel(int x, int y, Pixel p);
	void outputPPM(char *filename);
};

ColorImage::ColorImage()
{
	pPixel = 0; //指定pointer為0 即初始化Null
}

ColorImage::~ColorImage()
{
	if (pPixel) delete[] pPixel;
	pPixel = 0; //指定pointer為0 即初始化Null
}

void ColorImage::init(int xSize, int ySize)
{
	Pixel p = {0,0,0};
	xRes = xSize;
	yRes = ySize;
	pPixel = new Pixel[xSize*ySize];
	clear(p);
}

void ColorImage::clear(Pixel background) //background color
{
	int i;

	if (! pPixel) return;
	for (i=0; i<xRes*yRes; i++) pPixel[i] = background;
}

Pixel ColorImage::readPixel(int x, int y)
{
	assert(pPixel); // die if image not initialized
	return pPixel[x + y*yRes];
}

void ColorImage::writePixel(int x, int y, Pixel p)
{
	assert(pPixel); // die if image not initialized
	pPixel[x + y*xRes] = p;
}

void ColorImage::outputPPM(char *filename)
{
    FILE *outFile = fopen(filename, "wb");

	assert(outFile); // die if file can't be opened

	fprintf(outFile, "P6 %d %d 255\n", xRes, yRes);
	fwrite(pPixel, 1, 3*xRes*yRes, outFile );

	fclose(outFile);
}
//------------------------create ppm image functon End-----------------------------------------



#endif