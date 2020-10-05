#ifndef __GRAB_H__
#define __GRAB_H__

#include <iostream>
#include <fstream>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

#include "webcam.h"
#include "raylib.h"
//#include <MathCore.cpp>

#define XRES (800)
#define YRES (600)

Color toGreyScale(Color);
void makeScreenshot(const RGBImage*, const char*);

struct Frame
{
    size_t w;
    size_t h;
    Color *data;

    void draw()
    {
     for(size_t i=0;i<h;i++)
        {
              const int index = i*w;
            for(size_t j=0;j<w;j++)
            {
                const register int jj = index+j;
                DrawPixel(j,i,data[jj]);
            }
        }
    }


    Frame(const RGBImage &image)
    {
        data = new Color[image.height*image.width];
        w = image.width;
        h = image.height;
        for(size_t i=0;i<image.height;i++)
        {
              const int index = i*w;
            for(size_t j=0;j<image.width;j++)
            {
                const register int jj = index+j;
                size_t k = 3*(i*image.width+j);
                data[jj].r = image.data[k];
                data[jj].g = image.data[k+1];
                data[jj].b = image.data[k+2];
                data[jj].a = 255;
            }
        }
     }
    ~Frame()
    {
        delete[] data;
     }
};


#endif
