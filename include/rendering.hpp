#pragma once

void rgb_split(unsigned long c, float *r, float *g, float *b);
unsigned long rgb_join(float r, float g, float b);
void ppm_set(unsigned char *buf, int x, int y, unsigned long color);
unsigned long ppm_get(unsigned char *buf, int x, int y);
void ppm_write(const unsigned char *buf, FILE *f);
void pgm_write(const unsigned char *buf, std::FILE *f);
void frame(const real *var);
