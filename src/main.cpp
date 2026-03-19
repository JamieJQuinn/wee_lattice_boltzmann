#include <iostream>
#include <chrono>
#include <cmath>

#include "parameters.hpp"
#include "rendering.hpp"

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

inline int idx (const int i, const int j) {return i*NY + j;}
inline int idx (const int i, const int j, const int k) {return idx(i,j)*NUM_SPEEDS + k;}

inline real calc_f_eq_kernel(const real rho, const real u, const real v, int k, const real cs2) {
  real c_dot_u = C_VECS[k][0] * u + C_VECS[k][1] * v;
  real u2 = u*u + v*v;
  return WEIGHTS[k] * rho * (1.0 + c_dot_u/cs2 + 0.5*c_dot_u*c_dot_u/(cs2*cs2) - 0.5*u2/cs2);
}

void calc_f_eq(real* f_eq, const real* rho, const real* u, const real* v, const real cs2) {
#pragma omp parallel for collapse(3)
  for(int i=0; i<NX; ++i) {
    for (int j=0; j<NY; ++j) {
      for(int k=0; k<NUM_SPEEDS; ++k) {
        f_eq[idx(i,j,k)] = calc_f_eq_kernel(rho[idx(i,j)], u[idx(i,j)], v[idx(i,j)], k, cs2);
      }
    }
  }
}

void stream(real* f_new, const real* f_old) {
#pragma omp parallel for collapse(3)
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

void collide(real* f, const real* f_eq, const real inv_tau) {
#pragma omp parallel for collapse(3)
  for(int i=0; i<NX; ++i) {
    for (int j=0; j<NY; ++j) {
      for(int k=0; k<NUM_SPEEDS; ++k) {
        f[idx(i,j,k)] = (1.0 - inv_tau)*f[idx(i, j, k)] + inv_tau*f_eq[idx(i, j, k)];
      }
    }
  }
}

void calc_bounce_back(real* f_bnd, const real* f, const bool *obstacle) {
#pragma omp parallel for collapse(2)
  for(int i=0; i<NX; ++i) {
    for (int j=0; j<NY; ++j) {
      if (obstacle[idx(i,j)]) {
        for(int k=0; k<NUM_SPEEDS; ++k) {
          f_bnd[idx(i,j,k)] = f[idx(i,j,k)];
        }
        for(int i_bounce=0; i_bounce<4; ++i_bounce) {
          std::swap(f_bnd[idx(i,j,BOUNCE_BACKS[i_bounce][0])], f_bnd[idx(i,j,BOUNCE_BACKS[i_bounce][1])]);
        }
      }
    }
  }
}

void apply_bounce_back(real* f, const real* f_bnd, const bool *obstacle) {
#pragma omp parallel for collapse(2)
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

void update_macro_vars(real* rho, real* u, real* v, const real* f) {
#pragma omp parallel for collapse(2)
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
  if(tau < 0.5) {
    std::cout << "tau (" << tau << ") < 0.5 and is unstable. Stopping." << std::endl;
    exit(-1);
  }

  std::cout << "Re: " << Re << std::endl;
  std::cout << "u_pipe: " << u_pipe << std::endl;
  std::cout << "visc: " << visc << std::endl;
  const real Re_grid = (u_pipe/u0)/visc;
  std::cout << "Initial Re_grid: " << Re_grid << std::endl;
  std::cout << "tau: " << tau << std::endl;

  std::cout << "NX: " << NX << std::endl;
  std::cout << "NY: " << NY << std::endl;
  std::cout << "dx: " << dx << std::endl;
  std::cout << "dt: " << dt << std::endl;
  std::cout << "total_time: " << total_time << std::endl;
  const int total_steps = int(total_time/dt);
  std::cout << "n_steps: " << total_steps << std::endl;

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
      rho[idx(i,j)] = rho0;
      u[idx(i,j)] = u_pipe/u0;
      v[idx(i,j)] = 0.0/u0;
    }
  }
  calc_f_eq(f, rho, u, v, cs2);

  // Add some jiggle
  for(int i=0; i<NX; ++i) {
    for(int j=0; j<NY; ++j) {
      for(int k=0; k<NUM_SPEEDS; ++k) {
        f[idx(i,j,k)] += 0.01*(cos(M_PI*13.13*i) + sin(M_PI*19.67*j) + cos(M_PI*M_PI*k));
      }
    }
  }

  // DEFINE OBSTACLE

  int cx = NX/8; int cy = NY/2; int radius = NY/16;
  bool *cylinder = new bool[NX*NY];

  for(int i=0; i<NX; ++i) {
    for (int j=0; j<NY; ++j) {
      cylinder[idx(i,j)] = (pow(i-cx, 2) + pow(j-cy,2) < pow(radius, 2));
    }
  }

  // MAIN LOOP

  real t = 0.0;
  real t_until_next_frame = 0;
  real max = 0.0;
  real last_max = 0.0;
  unsigned int counter = 0;
  unsigned int dump_counter = 0;

  auto start = high_resolution_clock::now();
  while(t < total_time) {
    // apply_forces(...)
    update_macro_vars(rho, u, v, f);
    calc_f_eq(f_eq, rho, u, v, cs2);
    collide(f, f_eq, inv_tau);
    apply_bounce_back(f, f_bnd, cylinder);
    stream(temp, f);
    std::swap(f, temp);
    apply_f_bcs(f);
    calc_bounce_back(f_bnd, f, cylinder);

    // Apply pressure difference in x
    for (int j=0; j<NY; ++j) {
      for(int k=0; k<NUM_SPEEDS; ++k) {
        f[idx(0,j,k)] += calc_f_eq_kernel(rho[idx(0,j)] + rho_diff, u[idx(0,j)], v[idx(0,j)], k, cs2) - f_eq[idx(0,j,k)];
        f[idx(NX-1,j,k)] += calc_f_eq_kernel(rho[idx(NX-1,j)], u[idx(NX-1,j)], v[idx(NX-1,j)], k, cs2) - f_eq[idx(NX-1,j,k)];
      }
    }

    if (t > t_until_next_frame) {
      t_until_next_frame += dt_frame;
      for(int i=1; i<NX-1; ++i) {
        for(int j=1; j<NY-1; ++j) {
          if(cylinder[idx(i,j)]) {
            u[idx(i,j)] = 0.0;
            v[idx(i,j)] = 0.0;
          }
          real val = (v[idx(i+1,j)] - v[idx(i-1,j)])*NX - (u[idx(i,j+1)] - u[idx(i,j-1)])*NX;
          //real val = v[idx(i,j)];
          max = std::max(std::abs(val), max);
          real grey = ((val/last_max)*8 + 1.0)*0.5;
          ppm_set(ppm_buf, NX-i, j, rgb_join(grey, grey, grey));
        }
      }
      last_max = max;
      if(save_to_file) {
        std::string fname = format_counter(dump_counter) + ".ppm";
        std::FILE* fp = std::fopen(fname.c_str(), "w+");
        ppm_write(ppm_buf, fp, NX, NY);
        std::fclose(fp);
        dump_counter++;
      } else {
        ppm_write(ppm_buf, stdout, NX, NY);
      }

      real u_max = max_vel(u,v);
      if (u_max >= std::sqrt(2.0/3.0)) {
        std::cout << "u_max (= " << u_max << ") breached sqrt(2/3), increase resolution to compensate" << std::endl;
        exit(-1);
      }
    }

    t += dt;
    counter++;
}

  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - start).count();

  std::cout << "Time : " << duration << " ms" << std::endl;
  std::cout << "FPS : " << float(total_steps)/(duration/1000.0) << std::endl;
  std::cout << "time per step : " << duration/float(total_steps) << " ms" << std::endl;

  return 0;
}
