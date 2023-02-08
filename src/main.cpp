#include <cstring>
#include <iostream>
#include <chrono>
#include <cmath>
#include "parameters.hpp"
#include "rendering.hpp"

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

real WEIGHTS[NUM_SPEEDS] = {
  4.0 / 9.0,
  1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
  1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
int C_VECS[NUM_SPEEDS][2] = {
  {0, 0}, 
  {1, 0}, {0, 1}, {-1, 0}, {0, -1},
  {1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
int OPPOSITES[NUM_SPEEDS] = {
  0,
  3, 4, 1, 2,
  7, 8, 5, 6};
int BOUNCE_BACKS[4][2] = {
  {1, 3}, {2, 4}, {5, 7}, {6, 8}
};

inline int idx (const int i, const int j, const int k) {
  return idx(i,j)*NUM_SPEEDS + k;
}

real calc_f_eq_kernel(const real rho, const real u, const real v, int k, const real cs2) {
  real c_dot_u = C_VECS[k][0] * u + C_VECS[k][1] * v;
  real u2 = u*u + v*v;
  return WEIGHTS[k] * rho * (1.0 + c_dot_u/cs2 + 0.5*c_dot_u*c_dot_u/(cs2*cs2) - 0.5*u2/cs2);
}

void calc_f_eq(real *f_eq, const real *rho, const real *u, const real *v, const real cs2) {
  for(int i=0; i<NX; ++i) {
    for (int j=0; j<NY; ++j) {
      for(int k=0; k<NUM_SPEEDS; ++k) {
        f_eq[idx(i,j,k)] = calc_f_eq_kernel(rho[idx(i,j)], u[idx(i,j)], v[idx(i,j)], k, cs2);
      }
    }
  }
}

void stream(real *f_new, const real *f_old) {
  for(int i=1; i<NX-1; ++i) {
    for (int j=1; j<NY-1; ++j) {
      for(int k=0; k<NUM_SPEEDS; ++k) {
        int ip = i - C_VECS[k][0];
        int jp = j - C_VECS[k][1];
        f_new[idx(i,j,k)] = f_old[idx(ip, jp, k)];
      }
    }
  }
}

void collide(real *f, const real *f_eq, const real inv_tau) {
  for(int i=0; i<NX; ++i) {
    for (int j=0; j<NY; ++j) {
      for(int k=0; k<NUM_SPEEDS; ++k) {
        f[idx(i,j,k)] = (1.0 - inv_tau)*f[idx(i, j, k)] + inv_tau*f_eq[idx(i, j, k)];
      }
    }
  }
}

void apply_obstacle(real *f, const bool *obstacle) {
  for(int i=0; i<NX; ++i) {
    for (int j=0; j<NY; ++j) {
      if (obstacle[idx(i,j)]) {
        for(int i_bounce=0; i_bounce<4; ++i_bounce) {
          std::swap(f[idx(i,j,BOUNCE_BACKS[i_bounce][0])], f[idx(i,j,BOUNCE_BACKS[i_bounce][1])]);
        }
      }
    }
  }
}

void where(real *var, const bool *obstacle, const real val) {
  for(int i=0; i<NX; ++i) {
    for (int j=0; j<NY; ++j) {
      if (obstacle[idx(i,j)]) {
        var[idx(i,j)] = val;
      }
    }
  }
}

void update_macro_vars(real *rho, real *u, real *v, const real *f) {
  for(int i=0; i<NX; ++i) {
    for (int j=0; j<NY; ++j) {
      rho[idx(i,j)] = 0.0;
      u[idx(i,j)] = 0.0;
      v[idx(i,j)] = 0.0;
      for(int k=0; k<NUM_SPEEDS; ++k) {
        rho[idx(i,j)] += f[idx(i,j,k)];
        u[idx(i,j)] += C_VECS[k][0]*f[idx(i,j,k)];
        v[idx(i,j)] += C_VECS[k][1]*f[idx(i,j,k)];
      }
      u[idx(i,j)] /= rho[idx(i,j)];
      v[idx(i,j)] /= rho[idx(i,j)];
    }
  }
}

void apply_periodic_bcs(real *f) {
  for(int j=0; j<NY; ++j) {
    for(int k=0; k<NUM_SPEEDS; ++k) {
      f[idx(0,j,k)] = f[idx(NX-2,j,k)];
      f[idx(NX-1,j,k)] = f[idx(1,j,k)];
    }
  }

  for(int i=0; i<NX; ++i) {
    for(int k=0; k<NUM_SPEEDS; ++k) {
      f[idx(i,0,k)] = f[idx(i,NY-2,k)];
      f[idx(i,NY-1,k)] = f[idx(i,1,k)];
    }
  }
}

void apply_f_bcs(real *f) {
  apply_periodic_bcs(f);
}

std::string format_counter(int counter) {
  char buffer[5];
  std::snprintf(buffer, sizeof(buffer), "%04d", counter);
  return buffer;
}


real max_vel(const real *u, const real *v) {
  real max = pow(u[0], 2) + pow(v[0], 2);
  for(int i=0; i<NX; ++i) {
    for(int j=0; j<NY; ++j) {
      real val = pow(u[idx(i,j)], 2) + pow(v[idx(i,j)], 2);
      max = std::max(val, max);
    }
  }

  return std::sqrt(max);
}

int main() {
  const real LX = 1.0; // m
  const real dx = LX/NX; // m
  //const real dt = dx*dx; // diffusive scaling (s)
  const real dt = dx; // acoustic scaling (s)

  // nondimensionalisation
  const real u0 = dx/dt;
  const real rho0 = 1.0;

  const real Re = 100; // Reynold's Number
  const real u_pipe = 1.0;
  const real visc = LX*u_pipe/Re; // physical viscosity
  const real cs2 = 1.0/3.0; // nondmin sound speed
  //const real tau = visc/cs2 + 0.5;
  const real tau = 0.6;
  //const real tau = 0.55; // nondimensional
  const real inv_tau = 1.0 / tau;

  const real shear_visc = (cs2*(tau - 0.5));
  const real Re_grid = (u_pipe/u0)/shear_visc;
  //const real Re = u_pipe*LX/shear_visc;
  const real pressure_diff = 0.01;
  const real rho_diff = 0.01/cs2;

  const real total_time = 4000*dt; // s
  const real dt_frame = 100;

  const bool save_to_file = true;

  std::cerr << "Re: " << Re << std::endl;
  std::cerr << "Visc: " << shear_visc << std::endl;
  std::cerr << "Initial Re_grid: " << Re_grid << std::endl;
  //std::cerr << "lb_visc: " << lb_visc << std::endl;

  std::cerr << "NX: " << NX << std::endl;
  std::cerr << "NY: " << NY << std::endl;
  std::cerr << "dx: " << dx << std::endl;
  std::cerr << "dt: " << dt << std::endl;
  std::cerr << "tau: " << tau << std::endl;
  std::cerr << "total_time: " << total_time << std::endl;
  const int total_steps = int(total_time/dt);
  std::cerr << "n_steps: " << total_steps << std::endl;

  unsigned char *ppm_buf = new unsigned char[NX*NY*3];

  real *rho = new real[NX*NY];
  real *u = new real[NX*NY];
  real *v = new real[NX*NY];

  real *f_eq = new real[NX*NY*NUM_SPEEDS];
  real *f = new real[NX*NY*NUM_SPEEDS];
  real *temp = new real[NX*NY*NUM_SPEEDS];
  real *f_bnd = new real[NX*NY*NUM_SPEEDS];

  // INITIAL CONDITIONS

  for(int i=0; i<NX; ++i) {
    for(int j=0; j<NY; ++j) {
      rho[idx(i,j)] = 1.0/rho0;
      //u[idx(i,j)] = u_pipe/u0;
      u[idx(i,j)] = 2.0;
      v[idx(i,j)] = 0.0/u0;
    }
  }
  calc_f_eq(f, rho, u, v, cs2);

  // OBSTACLE

  int cx = NX/4; int cy = NY/2; int radius = NY/4;
  bool *cylinder = new bool[NX*NY];

  for(int i=0; i<NX; ++i) {
    for (int j=0; j<NY; ++j) {
      cylinder[idx(i,j)] = (pow(i-cx, 2) + pow(j-cy,2) <= pow(radius, 2));
    }
  }

  // MAIN LOOP

  //real max_vel = 0.0;
  real t = 0.0;
  real t_until_next_frame = 0;
  real max = 0.0;
  unsigned int counter = 0;
  unsigned int dump_counter = 0;

  auto start = high_resolution_clock::now();
  while(t < total_time) {
    // apply_forces(...)
    update_macro_vars(rho, u, v, f);
    calc_f_eq(f_eq, rho, u, v, cs2);
    collide(f, f_eq, inv_tau);
    stream(temp, f);
    std::swap(f, temp);
    apply_f_bcs(f);
    apply_obstacle(f, cylinder);

    // Apply pressure difference in x
    //for (int j=0; j<NY; ++j) {
      //for(int k=0; k<NUM_SPEEDS; ++k) {
        //f[idx(0,j,k)] += calc_f_eq_kernel(rho[idx(0,j)] + rho_diff, u[idx(0,j)], v[idx(0,j)], k) - f_eq[idx(0,j,k)];
        //f[idx(NX-1,j,k)] += calc_f_eq_kernel(rho[idx(NX-1,j)] + rho_diff, u[idx(NX-1,j)], v[idx(NX-1,j)], k) - f_eq[idx(NX-1,j,k)];
      //}
    //}

    // Update macro

    if (t > t_until_next_frame) {
      t_until_next_frame += dt_frame;
      for(int i=1; i<NX-1; ++i) {
        for(int j=1; j<NY-1; ++j) {
          real val = (v[idx(i+1,j)] - v[idx(i-1,j)])*NX - (u[idx(i,j+1)] - u[idx(i,j-1)])*NY;
          //real val = v[idx(i,j)];
          max = std::max(std::abs(val), max);
          real grey = ((val/max) + 1.0)/0.5;
          ppm_set(ppm_buf, NX-i, j, to_grey(grey));
        }
      }
      if(save_to_file) {
        std::string fname = format_counter(dump_counter) + ".ppm";
        std::FILE* fp = std::fopen(fname.c_str(), "w+");
        ppm_write(ppm_buf, fp);
        std::fclose(fp);
        dump_counter++;
      } else {
        ppm_write(ppm_buf, stdout);
      }

      //real Re_grid = max_vel(u,v)/(cs2*(tau - 0.5));
      //if (Re_grid > 10) {
        //std::cout << "Re_grid breached 10, increase resolution to compensate" << std::endl;
        //exit(-1);
      //}
      //for (int i=0; i<NX; ++i)  {
        //std::cout << int(rho[idx(i,NY/2)]) << " ";
      //}
      //std::cout << std::endl;
      //std::cout << std::endl;
    }

    t += dt;
    counter++;
}

  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - start).count();

  std::cerr << "Time : " << duration << " ms" << std::endl;
  std::cerr << "FPS : " << float(total_steps)/(duration/1000.0) << std::endl;
  std::cerr << "time per step : " << duration/float(total_steps) << " ms" << std::endl;

  return 0;
}
