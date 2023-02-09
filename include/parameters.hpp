#pragma once

const int NUM_SPEEDS = 9;
const int NX = 512;
const int NY = NX/8;

typedef float real;

inline int idx (const int i, const int j) {
  return i*NY + j;
}
