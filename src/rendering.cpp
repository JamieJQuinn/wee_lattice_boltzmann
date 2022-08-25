#include <cstring>
#include <cstdio>
#include <cmath>
#include "parameters.hpp"
#include "rendering.hpp"

void
rgb_split(unsigned long c, float *r, float *g, float *b)
{
    *r = ((c >> 16) / 255.0f);
    *g = (((c >> 8) & 0xff) / 255.0f);
    *b = ((c & 0xff) / 255.0f);
}

unsigned long
rgb_join(float r, float g, float b)
{
    unsigned long ir = std::roundf(r * 255.0f);
    unsigned long ig = std::roundf(g * 255.0f);
    unsigned long ib = std::roundf(b * 255.0f);
    return (ir << 16) | (ig << 8) | ib;
}

void
ppm_set(unsigned char *buf, int x, int y, unsigned long color)
{
    buf[idx(x,y) * 3 + 0] = color >> 16;
    buf[idx(x,y) * 3 + 1] = color >>  8;
    buf[idx(x,y) * 3 + 2] = color >>  0;
}

unsigned long
ppm_get(unsigned char *buf, int x, int y)
{
    unsigned long r = buf[idx(x,y) * 3 + 0];
    unsigned long g = buf[idx(x,y) * 3 + 1];
    unsigned long b = buf[idx(x,y) * 3 + 2];
    return (r << 16) | (g << 8) | b;
}

void
ppm_write(const unsigned char *buf, std::FILE *f)
{
    fprintf(f, "P6\n%d %d\n255\n", NY, NX);
    fwrite(buf, NY * 3, NX, f);
    fflush(f);
}

void
pgm_write(const unsigned char *buf, std::FILE *f)
{
    fprintf(f, "P5\n%d %d\n255\n", NX, NY);
    fwrite(buf, NX, NY, f);
    fflush(f);
}

unsigned char to_grey(const real val) {
  return rgb_join(val, val, val);
}


void
frame(const real *var)
{
    //static unsigned char buf[NX * NY * 3];
    //std::memset(buf, 0, sizeof(buf));
    //for(int i=0; i<NX; ++i) {
      //for(int j=0; j<NY; ++j) {
        //ppm_set(buf,i,j,to_grey(var[idx(i,j)]/0.57));
      //}
    //}
    //ppm_write(buf, stdout);
}
