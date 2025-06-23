static char help[] = "Solves a 1D transient heat problem.\n\n";
#include <iostream>
#include <cmath>
#include <vector>
#include <petsc.h>
#include "heat_solver.hpp"
#include "hdf5_tools.hpp"
#include "vtk_tools.hpp"
#include "Vec2Array.hpp"
#include <fstream>

int main(int argc,char **args)
{
    // PETSc Initialization
    PetscInitialize(&argc,&args,(char*)0,help);

    double start_time = MPI_Wtime();
    

    MPI_Comm        comm;
    PetscMPIInt     rank;

    comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);

    // Physical parameters
    double rho = 1.0, cp = 1.0, kappa = 1.0, L = 1.0; 

    // Meshing parameters
    int N = 50; 
    int M = 1;

    // Time stepping parameters
    double initial_time = 0.0, final_time = 1.0;
    double dt = (final_time - initial_time) / M;

    // Space meshing parameters
    double bc_left = 0.0, bc_right = L;
    double dx = (bc_right - bc_left) / N;

    // Stability factor
    PetscReal r1 = kappa * dt / (rho * cp * dx * dx);
    PetscReal r2 = dt / (rho * cp );
    
    // Stability check
    // if (r1 > 0.5 && rank == 0){
    //     throw std::runtime_error("warning: r  > 0.5, computation may be unstable\n");
    // }

    // Restart options
    bool is_restart = false;
    int restart_index = 0;
    double restart_time = 0.0;
    double restart_step = 1.0e-3;

    Vec            temp, temp_x, F;     // temperature solution, derivative of temperature, source vector
    PetscInt       i,col[3],rstart,rend,nlocal;
    PetscScalar    value[3];
    std::vector<double> x_coor(N+1, 0.0);


    // // Create vectors.  Note that we form 1 vector from scratch and
    // // then duplicate as needed. For this simple case let PETSc decide how
    // // many elements of the vector are stored on each processor. The second
    // // argument to VecSetSizes() below causes PETSc to decide.
    VecCreate(PETSC_COMM_WORLD, &temp);
    VecSetSizes(temp, PETSC_DECIDE, N);
    VecSetFromOptions(temp);
    VecDuplicate(temp, &temp_x);
    VecDuplicate(temp, &F);

    // // Identify the starting and ending mesh points on each
    // // processor for the interior part of the mesh.
    VecGetOwnershipRange(temp, &rstart, &rend);
    VecGetLocalSize(temp,&nlocal);

    // Set up the coordinates
    for (int ii = 0; ii < N+1; ii++){
        x_coor[ii] = ii * dx;
    }
    
    Manufactured_Solution * mn_sol = new Manufactured_Solution(rho, cp, kappa, r1, r2, 
             initial_time, final_time, dt, bc_left, bc_right, dx, N, M, x_coor);

    heat_solver * hsolver = new heat_solver(mn_sol);

    hdf5_tools *h5_tls = new hdf5_tools("SOL_TEMPERATURE.h5");

    Vec2Array *vec2arry = new Vec2Array(); 

    h5_tls->setup_hdf5();
    MPI_Barrier(PETSC_COMM_WORLD);

    // hsolver->explicitEuler(temp, F, h5_tls, vec2arry);
    hsolver->implicitEuler(temp, F, h5_tls, vec2arry);

    std::vector<std::vector<double>> temp_timeset;     // store the temperature arrays for all time steps
    h5_tls->read_h5("SOL_TEMPERATURE.h5", temp_timeset);   // read the h5 file that stores the temperature arrays
    MPI_Barrier(PETSC_COMM_WORLD);

    std::vector<std::vector<double>> exact_temp_timeset = mn_sol->get_exact_temp_timeset();  // obtain the exact temperature arrays for all time steps

    std::vector<double> temp_err(M+1); 
    if (rank == 0){
        double max_err = 0.0;
        for (int t = 0; t <= M; t++){
            for (int i = 1; i <= N; i++){
                double err = std::abs(temp_timeset[t][i-1] - exact_temp_timeset[t][i]);
                if (err > max_err) max_err = err;
            }
            temp_err[t] = max_err;
            max_err = 0.0;
        }
        std::cout << "errors:\n";
        for (int tt = 0; tt <= M; tt++){
            std::cout << "time t " << tt << ": " << temp_err[tt] << "\t";
            if (tt % 5 == 4) std::cout << std::endl;
        }
         std::cout << std::endl;

         // write the errors into binary file
         std::ofstream file("SOL_ERROR", std::ios::binary);
        if (!file.write(reinterpret_cast<const char*>(temp_err.data()), temp_err.size() * sizeof(double))){
            throw std::runtime_error("fail writing");
        }
        file.close();
    }
    MPI_Barrier(PETSC_COMM_WORLD);

    vtk_tools *vtk_tls = new vtk_tools();

    vtk_tls->write_vtk(temp_timeset, "SOL_TEMPERATURE", dx, dt, 1.0);
    // vtk_tls->write_vtk(exact_temp_timeset, "SOL_EXACT_TEMPERATURE", dx, dt, 1.0);
    MPI_Barrier(PETSC_COMM_WORLD);

    // Free heap-allocated objects
    delete mn_sol;
    delete hsolver;
    delete h5_tls;
    delete vec2arry;
    delete vtk_tls;

    // Free work space.  All PETSc objects should be destroyed when they
    // are no longer needed.
    VecDestroy(&temp); VecDestroy(&temp_x);
    VecDestroy(&F);

    double end_time = MPI_Wtime();
    PetscPrintf(PETSC_COMM_WORLD, "Total time: %.16f seconds\n", end_time - start_time);

    // Always call PetscFinalize() before exiting a program.  This routine
    //   - finalizes the PETSc libraries as well as MPI
    //   - provides summary and diagnostic information if certain runtime
    //     options are chosen (e.g., -log_view).
    PetscFinalize();
    return EXIT_SUCCESS;
}

// EOF
