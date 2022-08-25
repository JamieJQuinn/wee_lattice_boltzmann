#pragma once

const int NUM_SPEEDS = 9;
const int NX = 256;
const int NY = 64;

typedef float real;

inline int idx (const int i, const int j) {
  return i*NY + j;
}
