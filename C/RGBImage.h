

#ifndef __RGBIMAGE_H__
#define __RGBIMAGE_H__

struct RGBImage {
      unsigned char   *data; // RGB888 <=> RGB24
      size_t          width;
      size_t          height;
      size_t          size; // width * height * 3
};

#endif //__RGBIMAGE_H__
