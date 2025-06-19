#ifndef MANUFACTURED_SOLUTION
#define MANUFACTURED_SOLUTION

#include <petsc.h>
#include <iostream>
#include <cmath>
#include <vector>

class Manufactured_Solution
{
private:
public:
    // physical parameters
    double rho, cp, kappa, r1, r2, initial_time, final_time, dt, bc_left, bc_right, dx;
    int N, M;
    std::vector<double> &x_coor;
    // constructor
    Manufactured_Solution(double rho_val, double cp_val, double kappa_val, double r1_val, double r2_val, 
                          double initial_time_val, double final_time_val, double dt_val, 
                          double bc_left_val, double bc_right_val, double dx_val, 
                          int N_val, int M_val,
                          std::vector<double>  &x_coor_in)
        : rho(rho_val), cp(cp_val), kappa(kappa_val), r1(r1_val), r2(r2_val), 
          initial_time(initial_time_val), final_time(final_time_val), dt(dt_val),
          bc_left(bc_left_val), bc_right(bc_right_val), dx(dx_val),
          N(N_val), M(M_val), x_coor(x_coor_in) {}
    
    // destructor
    ~Manufactured_Solution();

    // heat supply per unit volume function: rho c (2t-3)(2-x)(3-x) - 2 kappa (t-2)(t-1)
    double f(double x, double t) const {
        // return rho * cp * (2*t - 3) * (2 - x) * (3 - x) - 2 * kappa * (t - 2) * (t - 1);
        return 0.1 * (rho * cp * cos(t) * (2 - x) * (3 - x) - 2 * kappa * sin(t));
    }

    // Dirichlet B.C. (t-2)(t-1)(2-x)(3-x) 
    double g(double x, double t) const {
        // return (t - 2) * (t - 1) * (2 - x) * (3 - x);
        return 0.1 * sin(t) * (2 - x) * (3 - x);
    }

    // Neumann B.C. kappa n_x (t-2)(t-1)(2x-5)
    double h(double x, double t, double nx) const {
        // return kappa * nx * (t - 2) * (t - 1) * (2*x - 5);
        return 0.1 * kappa * nx * sin(t) * (2*x - 5);
    }

    // temperature field. (t-2)(t-1)(2-x)(3-x) 
    double u(double x, double t) const {
        // return (t - 2) * (t - 1) * (2 - x) * (3 - x);
        return 0.1 * sin(t) * (2 - x) * (3 - x);
    }

    // get exact solution for temperature field
    std::vector<std::vector<double>> get_exact_temp_timeset(){

        std::vector<std::vector<double>> temp_timeset(M+1, std::vector<double>(N+1, 0.0));
        double x_val, t_val;
        for (int tt = 0; tt <= M; tt++){
            t_val = tt * dt;
            for (int ii = 0; ii <= N; ii++){
                x_val = ii*dx;
                temp_timeset[tt][ii] = u(x_val, t_val);
            }
        }
        return temp_timeset;
    }

};




#endif