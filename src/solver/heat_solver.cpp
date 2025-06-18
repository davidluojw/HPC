// heat_solver.cpp
#include "heat_solver.hpp"
#include <petsc.h>

heat_solver::heat_solver(Manufactured_Solution *mn_sol_in)
    : mn_sol(mn_sol_in) {}


void heat_solver::initialize(Vec &temp, Vec &F) {
    PetscInt rstart, rend;
    PetscMPIInt    rank,size;
    VecGetOwnershipRange(temp, &rstart, &rend);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    double t0 = 0.0;

    // set temp = u(x,0)
    for (PetscInt i = rstart; i < rend; i++) {
        // manufactured solution u at t = 0:  u(x,0)
        double u_val = mn_sol->u(mn_sol->x_coor[i+1], t0); // internal nodes start from x_coor = dx
        
        // set Vec values
        VecSetValues(temp, 1, &i, &u_val, INSERT_VALUES);
    }

    // set F = f(x,0)
    for (PetscInt i = rstart; i < rend; i++) {        
        // manufactured solution f at t = 0: f(x,0)
        double f_val = mn_sol->r2 * mn_sol->f(mn_sol->x_coor[i+1], t0); // internal nodes start from x_coor = dx
        PetscPrintf(PETSC_COMM_SELF,"from rank %d , rstart = %d, rend = %d, i = %d, f_val = %f\n",rank, rstart, rend, i, f_val);
        if (i == 0){
            f_val += mn_sol->r1 * mn_sol->g(mn_sol->x_coor[0], t0);
            PetscPrintf(PETSC_COMM_SELF,"from rank %d , i = %d, f_val = %f\n",rank, i, f_val);
        }
        else if (i == mn_sol->N - 1){
            f_val += mn_sol->r1 * 2 * mn_sol->dx / mn_sol->kappa * mn_sol->h(mn_sol->x_coor[mn_sol->N], t0, 1);
            PetscPrintf(PETSC_COMM_SELF,"from rank %d, i = %d, f_val = %f \n",rank, i, f_val);
        }
        
        // set Vec values
        VecSetValues(F, 1, &i, &f_val, INSERT_VALUES);
    }

    // Assembly Vec
    VecAssemblyBegin(temp);
    VecAssemblyEnd(temp);
    VecAssemblyBegin(F);
    VecAssemblyEnd(F);

    
}

// void HeatSolver::timeStepLoop() {
//     int steps = T_ / dt_;
//     for (int step = 0; step <= steps; ++step) {
//         stepExplicit();
//         VecCopy(u_new_, u_);
//     }
// }

// void HeatSolver::stepExplicit() {
//     PetscScalar *u_array, *u_new_array;
//     PetscInt nlocal;
//     VecGetLocalSize(u_, &nlocal);
//     VecGetArray(u_, &u_array);
//     VecGetArray(u_new_, &u_new_array);

//     double dx = 1.0 / (nx_ - 1);
//     double dy = 1.0 / (ny_ - 1);
//     double coef = dt_ * kappa_ / (rho_ * c_);

//     for (int j = 0; j < ny_; ++j) {
//         for (int i = 0; i < nx_; ++i) {
//             int idx = j * nx_ + i;
//             double uij = u_array[idx];

//             // neighbors
//             double uL = (i > 0)     ? u_array[idx - 1] : 0.0;
//             double uR = (i < nx_ - 1) ? u_array[idx + 1] : 0.0;
//             double uD = (j > 0)     ? u_array[idx - nx_] : 0.0;
//             double uU = (j < ny_ - 1) ? u_array[idx + nx_] : 0.0;

//             double d2udx2 = (uL - 2*uij + uR) / (dx*dx);
//             double d2udy2 = (uD - 2*uij + uU) / (dy*dy);

//             u_new_array[idx] = uij + coef * (d2udx2 + d2udy2);
//         }
//     }

//     VecRestoreArray(u_, &u_array);
//     VecRestoreArray(u_new_, &u_new_array);
// }
