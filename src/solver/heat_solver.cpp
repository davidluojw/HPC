// heat_solver.cpp
#include "heat_solver.hpp"
#include "hdf5_tools.hpp"
#include "Vec2Array.hpp"
#include <petsc.h>

heat_solver::heat_solver(Manufactured_Solution *mn_sol_in)
    : mn_sol(mn_sol_in) {}


void heat_solver::initialize(Vec &temp, Vec &F) {
    PetscInt rstart, rend;
    PetscMPIInt    rank,size;
    VecGetOwnershipRange(temp, &rstart, &rend);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    double t0 = mn_sol->initial_time;

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

    VecView(temp, PETSC_VIEWER_STDOUT_WORLD);
    VecView(F, PETSC_VIEWER_STDOUT_WORLD);

}

void heat_solver::update_temperature_expliciteuler(Vec &temp, Vec &F, Mat &A){
    Vec temp_np1;
    VecDuplicate(temp, &temp_np1);

    // temp_np1 = A*temp
    MatMult(A, temp, temp_np1);

    // temp_np1 = temp_np1 + F
    VecAXPY(temp_np1, 1.0, F);

    // Update temp
    VecSwap(temp, temp_np1);

    // release
    VecDestroy(&temp_np1);

}

void heat_solver::update_sourcevec(Vec &F, int j){
    PetscInt rstart, rend;
    PetscMPIInt    rank,size;
    VecGetOwnershipRange(F, &rstart, &rend);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    double tt = j * mn_sol->dt;

    // set F = f(x,0)
    for (PetscInt i = rstart; i < rend; i++) {        
        // manufactured solution f at t = 0: f(x,0)
        double f_val = mn_sol->r2 * mn_sol->f(mn_sol->x_coor[i+1], tt); // internal nodes start from x_coor = dx
        PetscPrintf(PETSC_COMM_SELF,"from rank %d , rstart = %d, rend = %d, i = %d, f_val = %f\n",rank, rstart, rend, i, f_val);
        if (i == 0){
            f_val += mn_sol->r1 * mn_sol->g(mn_sol->x_coor[0], tt);
            PetscPrintf(PETSC_COMM_SELF,"from rank %d , i = %d, f_val = %f\n",rank, i, f_val);
        }
        else if (i == mn_sol->N - 1){
            f_val += mn_sol->r1 * 2 * mn_sol->dx / mn_sol->kappa * mn_sol->h(mn_sol->x_coor[mn_sol->N], tt, 1);
            PetscPrintf(PETSC_COMM_SELF,"from rank %d, i = %d, f_val = %f \n",rank, i, f_val);
        }
        
        // set Vec values
        VecSetValues(F, 1, &i, &f_val, INSERT_VALUES);
    }

    // Assembly Vec
    VecAssemblyBegin(F);
    VecAssemblyEnd(F);
}

void heat_solver::time_loop(Vec &temp, Vec &F, Mat &A, hdf5_tools * const & h5_tls, Vec2Array * const & vec2arry) {
    PetscInt rstart, rend;
    PetscMPIInt    rank,size;
    VecGetOwnershipRange(F, &rstart, &rend);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    for (int j = 1; j <= 2; j++) {
        update_temperature_expliciteuler(temp, F, A);
        update_sourcevec(F, j);

        VecView(temp, PETSC_VIEWER_STDOUT_WORLD);
        VecView(F, PETSC_VIEWER_STDOUT_WORLD);

        std::vector<double> step_temp = vec2arry->get_vector_array(temp);
        if (rank == 0){
            std::cout << "step_data: time_step = " << j << "\n";
            for (int ii = 0; ii < mn_sol->N; ++ii){
                std::cout << step_temp[ii] << "\t";
            }
            std::cout << std::endl;
        }

        h5_tls->write_hdf5(j, 0, step_temp);

    }

}
