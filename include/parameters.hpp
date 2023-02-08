#pragma once

const int NUM_SPEEDS = 9;
const int NX = 400;
const int NY = NX/4;

typedef float real;

inline int idx (const int i, const int j) {
  return i*NY + j;
}
