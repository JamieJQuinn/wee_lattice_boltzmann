#include <cstring>
#include <iostream>
#include <chrono>
#include <cmath>
#include "parameters.hpp"
#include "rendering.hpp"

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

real WEIGHTS[NUM_SPEEDS] = {4.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 9.0, 1.0 / 36.0};
int C_VECS[NUM_SPEEDS][2] = {{0, 0}, {1, 0}, {1, -1}, {0, -1}, {-1, -1}, {-1, 0}, {-1, 1}, {0, 1}, {1, 1}};
int OPPOSITES[NUM_SPEEDS] = {0, 5, 6, 7, 8, 1, 2, 3, 4};

inline int idx (const int i, const int j, const int k) {
  return idx(i,j)*NUM_SPEEDS + k;
}

void calc_f_eq(real *f_eq, const real *rho, const real *u, const real *v) {
  for(int i=0; i<NX; ++i) {
    for (int j=0; j<NY; ++j) {
      for(int k=0; k<NUM_SPEEDS; ++k) {
        real c_dot_u = C_VECS[k][0] * u[idx(i,j)] + C_VECS[k][1] * v[idx(i,j)];
        real u2 = u[idx(i,j)]*u[idx(i,j)] + v[idx(i,j)]*v[idx(i,j)];
        f_eq[idx(i,j,k)] = WEIGHTS[k] * rho[idx(i,j)] * (1.0 + 3.0*c_dot_u + 4.5*c_dot_u*c_dot_u - 1.5*u2);
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

void calc_obstacle(real *f_bnd, const real *f, const bool *obstacle) {
  for(int i=0; i<NX; ++i) {
    for (int j=0; j<NY; ++j) {
      if (obstacle[idx(i,j)]) {
        for(int k=0; k<NUM_SPEEDS; ++k) {
          f_bnd[idx(i,j,k)] = f[idx(i,j,OPPOSITES[k])];
        }
      }
    }
  }
}

void apply_obstacle(real *f, const real *f_bnd, const bool *obstacle) {
  for(int i=0; i<NX; ++i) {
    for (int j=0; j<NY; ++j) {
      if (obstacle[idx(i,j)]) {
        for(int k=0; k<NUM_SPEEDS; ++k) {
          f[idx(i,j,k)] = f_bnd[idx(i,j,k)];
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

int main() {
  //const real Re = 1000; // Reynold's Number
  const real dx = 1.0/NX;
  const real dt = 2*dx*dx;
  const real rho0 = 100.0; // average density

  //const real cs2 = 1.0/3.0; // Sound speed

  //const real lb_visc = dt/(dx*dx)/Re;
  //const real tau = lb_visc/cs2 + 0.5;
  const real tau = 0.6;
  const real inv_tau = 1.0 / tau;

  //const int n_frames = 100;
  //const real total_time = n_frames * dt;
  const real total_time = 3000*dt;
  const real dt_frame = dt;

  //std::cerr << "Re: " << Re << std::endl;
  std::cerr << "NX: " << NX << std::endl;
  //std::cerr << "dx: " << dx << std::endl;
  //std::cerr << "dt: " << dt << std::endl;
  //std::cerr << "lb_visc: " << lb_visc << std::endl;
  std::cerr << "tau: " << tau << std::endl;
  std::cerr << "total_time: " << total_time << std::endl;
  const int total_steps = int(total_time/dt);
  std::cerr << "n_steps: " << total_steps << std::endl;

  //const double time_per_step = 0.111; // s per 512^2 domain on mobile RTX 3060
  //const double time_to_completion = NX*NY/(512.0*512.0)*time_per_step*total_steps;
  //std::cerr << "time_to_completion (s): " << time_to_completion << std::endl;

  unsigned char *ppm_buf = new unsigned char[NX*NY*3];

  real *rho = new real[NX*NY];
  real *u = new real[NX*NY];
  real *v = new real[NX*NY];

  real *f_eq = new real[NX*NY*NUM_SPEEDS];
  real *f = new real[NX*NY*NUM_SPEEDS];
  real *temp = new real[NX*NY*NUM_SPEEDS];
  real *f_bnd = new real[NX*NY*NUM_SPEEDS];

  for(int i=0; i<NX; ++i) {
    for(int j=0; j<NY; ++j) {
      for(int k=0; k<NUM_SPEEDS; ++k) {
        f[idx(i,j,k)] = 1.0 + 0.01*sin((i*j)*dx*100);
      }
      f[idx(i,j,1)] += 2. * (1.+0.2*cos(2.*M_PI*float(i)/NX*4));
    }
  }

  update_macro_vars(rho, u, v, f);

  for(int i=0; i<NX; ++i) {
    for(int j=0; j<NY; ++j) {
      for(int k=0; k<NUM_SPEEDS; ++k) {
        f[idx(i,j,k)] *= rho0/rho[idx(i,j)];
      }
    }
  }

  int cx = NX/4; int cy = NY/2; int radius = NY/4;
  bool *cylinder = new bool[NX*NY];

  for(int i=0; i<NX; ++i) {
    for (int j=0; j<NY; ++j) {
      cylinder[idx(i,j)] = (pow(i-cx, 2) + pow(j-cy,2) <= pow(radius, 2));
    }
  }

  auto start = high_resolution_clock::now();

  real t = 0.000001;
  real t_until_next_frame = 0;
  real max = 0.0;
  {
    while(t < total_time) {
      apply_periodic_bcs(f);
      stream(temp, f);
      std::swap(f, temp);
      apply_periodic_bcs(f);
      calc_obstacle(f_bnd, f, cylinder);
      update_macro_vars(rho, u, v, f);
      calc_f_eq(f_eq, rho, u, v);
      collide(f, f_eq, inv_tau);
      apply_obstacle(f, f_bnd, cylinder);
      where(u, cylinder, 0.0);
      where(v, cylinder, 0.0);

      //if (t > t_until_next_frame) {
        //t_until_next_frame += dt_frame;
        //for(int i=1; i<NX-1; ++i) {
          //for(int j=1; j<NY-1; ++j) {
            ////real val = (v[idx(i+1,j)] - v[idx(i-1,j)])*NX - (u[idx(i,j+1)] - u[idx(i,j-1)])*NX;
            //real val = v[idx(i,j)];
            //max = std::max(std::abs(val), max);
            //real grey = ((val/max) + 1.0)*0.5;
            //ppm_set(ppm_buf, i, j, rgb_join(grey,grey,grey));
          //}
        //}
        //ppm_write(ppm_buf, stdout);
        ////for (int i=0; i<NX; ++i)  {
          ////std::cout << int(rho[idx(i,NY/2)]) << " ";
        ////}
        ////std::cout << std::endl;
        ////std::cout << std::endl;
      //}

      t += dt;
    }
  }

  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - start).count();

  std::cerr << "Time : " << duration << " ms" << std::endl;
  std::cerr << "FPS : " << float(total_steps)/(duration/1000.0) << std::endl;
  std::cerr << "time per step : " << duration/float(total_steps) << " ms" << std::endl;

  return 0;
}
