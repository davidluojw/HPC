static char help[] = "Solves a 1D transient heat problem.\n\n";
#include <iostream>
#include <cmath>
#include <vector>
#include <petsc.h>
#include "heat_solver.hpp"
#include "hdf5_tools.hpp"
#include "vtk_tools.hpp"
#include "Vec2Array.hpp"

int main(int argc,char **args)
{
    // PETSc Initialization
    PetscInitialize(&argc,&args,(char*)0,help);

    MPI_Comm        comm;
    PetscMPIInt     rank;

    comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);

    // Physical parameters
    double rho = 1.0, cp = 1.0, kappa = 1.0, L = 1.0; 

    // Meshing parameters
    int N = 10; 
    int M = 2000;

    // Time stepping parameters
    double initial_time = 0.0, final_time = 5.0;
    double dt = (final_time - initial_time) / M;

    // Space meshing parameters
    double bc_left = 0.0, bc_right = L;
    double dx = (bc_right - bc_left) / N;

    // Stability factor
    PetscReal r1 = kappa * dt / (rho * cp * dx * dx);
    PetscReal r2 = dt / (rho * cp );
    
    // Stability check
    if (r1 > 0.5 && rank == 0){
        throw std::runtime_error("warning: r  > 0.5, computation may be unstable\n");
    }

    // Restart options
    bool is_restart = false;
    int restart_index = 0;
    double restart_time = 0.0;
    double restart_step = 1.0e-3;

    Vec            temp, temp_x, F;     // temperature solution, derivative of temperature, source vector
    Mat            A;                   // linear system matrix
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
    if (rank == 0){
        std::cout << "x_coor: \n";
        for (int ii = 0; ii < N+1; ++ii){
            std::cout << x_coor[ii] << "\t";
        }
        std::cout << std::endl;
    }
    
    // // Create matrix.
    // // We pass in nlocal as the "local" size of the matrix to force it
    // // to have the same parallel layout as the vector created above.
    MatCreate(PETSC_COMM_WORLD,&A);
    MatSetSizes(A, nlocal, nlocal, N, N);
    MatSetType(A, MATMPIAIJ); 
    MatSetFromOptions(A);
    MatSetUp(A);

    // // Assemble matrix.
    // // The linear system is distributed across the processors by
    // // chunks of contiguous rows, which correspond to contiguous
    // // sections of the mesh on which the problem is discretized.
    if (rstart == 0) 
    {
        rstart = 1;
        i      = 0; col[0] = 0; col[1] = 1; value[0] = 1 - 2.0*r1; value[1] = r1;
        MatSetValues(A, 1, &i, 2, col, value, INSERT_VALUES);
    }

    if (rend == N) 
    {
        rend = N-1;
        i    = N-1; col[0] = N-2; col[1] = N-1; value[0] = 2.0*r1; value[1] = 1 - 2.0*r1;
        MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);
    }

    // Set entries corresponding to the mesh interior
    value[0] = r1; value[1] = 1 - 2.0*r1; value[2] = r1;
    for (i=rstart; i<rend; i++) 
    {
        col[0] = i-1; col[1] = i; col[2] = i+1;
        MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);
    }

    // Assemble the matrix
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

    MatView(A,PETSC_VIEWER_STDOUT_WORLD);

    Manufactured_Solution * mn_sol = new Manufactured_Solution(rho, cp, kappa, r1, r2, 
             initial_time, final_time, dt, bc_left, bc_right, dx, N, M, x_coor);

    heat_solver * hsolver = new heat_solver(mn_sol);

    hdf5_tools *h5_tls = new hdf5_tools("SOL_TEMPERATURE.h5");

    Vec2Array *vec2arry = new Vec2Array(); 

    hsolver->initialize(temp, F); // initialize the temperature vector and source vector

    std::vector<double> init_temp = vec2arry->get_vector_array(temp);  // get the array of Vec temp
    if (rank == 0){
        std::cout << "init_temp: \n";
        for (int ii = 0; ii < N; ++ii){
            std::cout << init_temp[ii] << "\t";
        }
        std::cout << std::endl;
    }

    h5_tls->setup_hdf5();                    // create a h5 file
    h5_tls->write_hdf5(0, 0, init_temp);     // write the array into the h5 file

    hsolver->time_loop(temp, F, A, h5_tls, vec2arry);   // loop over the time step, meanwhile write the h5 file for each time step

    std::vector<std::vector<double>> temp_timesets;     // store the temperature arrays for all time steps
    h5_tls->read_h5("SOL_TEMPERATURE.h5", temp_timesets);   // read the h5 file that stores the temperature arrays
    if (rank == 0){
        std::cout << "temp_timesets: \n";
        for (int tt = 0; tt <= 2; tt++){
            std::cout << "time t " << tt << ": \t";
            for (int ii = 0; ii < N; ++ii){
                std::cout << temp_timesets[tt][ii] << "\t";
            }
            std::cout << std::endl;
        }
    }

    std::vector<std::vector<double>> exact_temp_timeset = mn_sol->get_exact_temp_timeset();  // obtain the exact temperature arrays for all time steps
    if (rank == 0){
        std::cout << "exact_temp_timesets: \n";
        for (int tt = 0; tt <= 2; tt++){
            std::cout << "time t " << tt << ": \t";
            for (int ii = 1; ii <= N; ++ii){
                std::cout << exact_temp_timeset[tt][ii] << "\t";
            }
            std::cout << std::endl;
        }
    }

    vtk_tools *vtk_tls = new vtk_tools();

    vtk_tls->write_vtk(temp_timesets, "SOL_TEMPERATURE", dx, dt, 1.0);
    vtk_tls->write_vtk(exact_temp_timeset, "SOL_EXACT_TEMPERATURE", dx, dt, 1.0);


    // // Set exact solution;
    // PetscScalar *u_array;
    // PetscScalar *b_array;
    // VecGetArray(u, &u_array);
    // VecGetArray(b, &b_array);
    // VecGetOwnershipRange(x, &rstart, &rend);
    // PetscReal pi = PETSC_PI;
    // for (i = rstart; i < rend; ++i) 
    // {
    //     PetscReal xi = (i + 1) * h;
    //     u_array[i - rstart] = sin(pi*xi);
    //     b_array[i - rstart] = h*h*pi*pi*sin(pi*xi);

    //     //PetscPrintf(PETSC_COMM_WORLD, "i = %d, xi = %e \n", i, xi);
    // }
    // VecRestoreArray(u, &u_array);
    // VecRestoreArray(b, &b_array);

    // //VecView(u,PETSC_VIEWER_STDOUT_WORLD);

    // // Create the linear solver and set various options
    // KSPCreate(PETSC_COMM_WORLD, &ksp);

    // // Set operators. Here the matrix that defines the linear system
    // // also serves as the preconditioning matrix.
    // KSPSetOperators(ksp, A, A);

    // // Set linear solver defaults for this problem (optional).
    // // - By extracting the KSP and PC contexts from the KSP context,
    // //   we can then directly call any KSP and PC routines to set
    // //   various options.
    // // - The following four statements are optional; all of these
    // //   parameters could alternatively be specified at runtime via
    // //   KSPSetFromOptions();
    // KSPGetPC(ksp, &pc);
    // PCSetType(pc, PCJACOBI);
    // KSPSetTolerances(ksp, 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

    // // Set runtime options, e.g.,
    // //     -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    // // These options will override those specified above as long as
    // // KSPSetFromOptions() is called _after_ any other customization routines.
    // KSPSetFromOptions(ksp);

    // // Solve the linear system
    // KSPSolve(ksp,b,x);

    // // Check solution and clean up
    // // Check the error
    // VecAXPY(x, -1.0, u);
    // VecNorm(x, NORM_2, &norm);
    // norm = norm * PetscSqrtReal(h);
    // KSPGetIterationNumber(ksp, &its);

    // PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g, Iterations %D\n",(double)norm,its);

    // // Free work space.  All PETSc objects should be destroyed when they
    // // are no longer needed.
    // VecDestroy(&x); VecDestroy(&u);
    // VecDestroy(&b); MatDestroy(&A);
    // KSPDestroy(&ksp);

    // // Always call PetscFinalize() before exiting a program.  This routine
    // //   - finalizes the PETSc libraries as well as MPI
    // //   - provides summary and diagnostic information if certain runtime
    // //     options are chosen (e.g., -log_view).
    // PetscFinalize();
    return 0;
}

// EOF
