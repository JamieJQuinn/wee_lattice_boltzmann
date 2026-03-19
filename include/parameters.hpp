#pragma once

typedef float real;

const int NY = 512;
const int NX = 1.777778*NY;

const real total_time = 20.; // s
const real dt_frame = total_time/500.;

const real u_pipe = 0.2; // initial velocity in pipe
const real pressure_diff = 0.00; // pressure difference applied across pipe

const bool save_to_file = true;

const int NUM_SPEEDS = 9;
const real WEIGHTS[NUM_SPEEDS] = {
  4.0 / 9.0,
  1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
  1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
const int C_VECS[NUM_SPEEDS][2] = {
  {0, 0}, 
  {1, 0}, {0, 1}, {-1, 0}, {0, -1},
  {1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
const int BOUNCE_BACKS[4][2] = {
  {1, 3}, {2, 4}, {5, 7}, {6, 8}
};

const real LX = 1.0; // m
const real dx = LX/NX; // m
//const real dt = dx*dx; // diffusive scaling (s)
const real dt = dx; // acoustic scaling (s)

// nondimensionalisation
const real u0 = dx/dt;
const real rho0 = 1.0;

const real Re = 10; // Reynold's Number
const real visc = LX*u_pipe/Re; // physical viscosity
const real cs2 = 1.0/3.0; // nondmin sound speed
//const real tau = visc/cs2 + 0.5;
const real tau = 0.55;
const real inv_tau = 1.0 / tau;

const real rho_diff = pressure_diff/cs2;
