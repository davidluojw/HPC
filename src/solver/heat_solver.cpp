// heat_solver.cpp
#include "heat_solver.hpp"
#include "hdf5_tools.hpp"
#include "Vec2Array.hpp"
#include <petsc.h>

heat_solver::heat_solver(Manufactured_Solution *mn_sol_in)
    : mn_sol(mn_sol_in) {}


void heat_solver::initialize(Vec &temp, Vec &F, hdf5_tools * const & h5_tls, Vec2Array * const & vec2arry) {
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
        // PetscPrintf(PETSC_COMM_SELF,"from rank %d , rstart = %d, rend = %d, i = %d, f_val = %f\n",rank, rstart, rend, i, f_val);
        if (i == 0){
            f_val += mn_sol->r1 * mn_sol->g(mn_sol->x_coor[0], t0);
            // PetscPrintf(PETSC_COMM_SELF,"from rank %d , i = %d, f_val = %f\n",rank, i, f_val);
        }
        else if (i == mn_sol->N - 1){
            f_val += mn_sol->r1 * 2 * mn_sol->dx / mn_sol->kappa * mn_sol->h(mn_sol->x_coor[mn_sol->N], t0, 1);
            // PetscPrintf(PETSC_COMM_SELF,"from rank %d, i = %d, f_val = %f \n",rank, i, f_val);
        }
        
        // set Vec values
        VecSetValues(F, 1, &i, &f_val, INSERT_VALUES);
    }

    // Assembly Vec
    VecAssemblyBegin(temp);
    VecAssemblyEnd(temp);
    VecAssemblyBegin(F);
    VecAssemblyEnd(F);

    // VecView(temp, PETSC_VIEWER_STDOUT_WORLD);
    // VecView(F, PETSC_VIEWER_STDOUT_WORLD);

    std::vector<double> init_temp = vec2arry->get_vector_array(temp);  // get the array of Vec temp
    if (rank == 0){
        // std::cout << "init_temp: time_step = 0 \n";
        // for (int ii = 0; ii < mn_sol->N; ++ii){
        //     std::cout << init_temp[ii] << "\t";
        // }
        // std::cout << std::endl;

        h5_tls->write_hdf5(0, 0, init_temp);
    }

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
        // PetscPrintf(PETSC_COMM_SELF,"from rank %d , rstart = %d, rend = %d, i = %d, f_val = %f\n",rank, rstart, rend, i, f_val);
        if (i == 0){
            f_val += mn_sol->r1 * mn_sol->g(mn_sol->x_coor[0], tt);
            // PetscPrintf(PETSC_COMM_SELF,"from rank %d , i = %d, f_val = %f\n",rank, i, f_val);
        }
        else if (i == mn_sol->N - 1){
            f_val += mn_sol->r1 * 2 * mn_sol->dx / mn_sol->kappa * mn_sol->h(mn_sol->x_coor[mn_sol->N], tt, 1);
            // PetscPrintf(PETSC_COMM_SELF,"from rank %d, i = %d, f_val = %f \n",rank, i, f_val);
        }
        
        // set Vec values
        VecSetValues(F, 1, &i, &f_val, INSERT_VALUES);
    }

    // Assembly Vec
    VecAssemblyBegin(F);
    VecAssemblyEnd(F);
}

void heat_solver::time_loop_expliciteuler(Vec &temp, Vec &F, Mat &A, hdf5_tools * const & h5_tls, Vec2Array * const & vec2arry) {
    PetscInt rstart, rend;
    PetscMPIInt    rank,size;
    VecGetOwnershipRange(F, &rstart, &rend);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    for (int j = 1; j <= mn_sol->M; j++) {
        update_temperature_expliciteuler(temp, F, A);
        update_sourcevec(F, j);

        // VecView(temp, PETSC_VIEWER_STDOUT_WORLD);
        // VecView(F, PETSC_VIEWER_STDOUT_WORLD);

        std::vector<double> step_temp = vec2arry->get_vector_array(temp);
        if (rank == 0){
            // std::cout << "step_temp: time_step = " << j << "\n";
            // for (int ii = 0; ii < mn_sol->N; ++ii){
            //     std::cout << step_temp[ii] << "\t";
            // }
            // std::cout << std::endl;

            h5_tls->write_hdf5(j, 0, step_temp);

        }

    }

}


void heat_solver::explicitEuler(Vec &temp, Vec &F, hdf5_tools * const & h5_tls, Vec2Array * const & vec2arry){
    int N = mn_sol->N;
    double r1 = mn_sol->r1;
    PetscInt       i,col[3],rstart,rend,nlocal;
    PetscScalar    value[3];

    // // Identify the starting and ending mesh points on each
    // // processor for the interior part of the mesh.
    VecGetOwnershipRange(temp, &rstart, &rend);
    VecGetLocalSize(temp,&nlocal);

    Mat            A;                   // linear system matrix

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

    initialize(temp, F, h5_tls, vec2arry); // initialize the temperature vector and source 

    time_loop_expliciteuler(temp, F, A, h5_tls, vec2arry);   // loop over the time step, meanwhile write the h5 file for each time step

     MatDestroy(&A);
}


void heat_solver::update_temperature_impliciteuler(Vec &temp, Vec &F, Mat &A){

    KSP ksp;
    PC pc;

    // F = temp + F
    VecAXPY(F, 1.0, temp);

     // Create the linear solver and set various options
    KSPCreate(PETSC_COMM_WORLD, &ksp);

    // Set operators. Here the matrix that defines the linear system
    // also serves as the preconditioning matrix.
    KSPSetOperators(ksp, A, A);

    // Set linear solver defaults for this problem (optional).
    // - By extracting the KSP and PC contexts from the KSP context,
    //   we can then directly call any KSP and PC routines to set
    //   various options.
    // - The following four statements are optional; all of these
    //   parameters could alternatively be specified at runtime via
    //   KSPSetFromOptions();
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCJACOBI);
    KSPSetTolerances(ksp, 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

    // Set runtime options, e.g.,
    //     -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    // These options will override those specified above as long as
    // KSPSetFromOptions() is called _after_ any other customization routines.
    KSPSetFromOptions(ksp);

    // Solve the linear system
    KSPSolve(ksp,F,temp);

    KSPDestroy(&ksp);

}

void heat_solver::time_loop_impliciteuler(Vec &temp, Vec &F, Mat &A, hdf5_tools * const & h5_tls, Vec2Array * const & vec2arry) {

    PetscInt rstart, rend;
    PetscMPIInt    rank,size;
    VecGetOwnershipRange(F, &rstart, &rend);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);


    for (int j = 1; j <= mn_sol->M; j++) {
        update_temperature_impliciteuler(temp, F, A);
        update_sourcevec(F, j+1);

        // VecView(temp, PETSC_VIEWER_STDOUT_WORLD);
        // VecView(F, PETSC_VIEWER_STDOUT_WORLD);

        std::vector<double> step_temp = vec2arry->get_vector_array(temp);
        if (rank == 0){
            // std::cout << "step_temp: time_step = " << j << "\n";
            // for (int ii = 0; ii < mn_sol->N; ++ii){
            //     std::cout << step_temp[ii] << "\t";
            // }
            // std::cout << std::endl;

            h5_tls->write_hdf5(j, 0, step_temp);

        }

    }

}

void heat_solver::implicitEuler(Vec &temp, Vec &F, hdf5_tools * const & h5_tls, Vec2Array * const & vec2arry){
    int N = mn_sol->N;
    double r1 = mn_sol->r1;
    PetscInt       i,col[3],rstart,rend,nlocal;
    PetscScalar    value[3];

    // // Identify the starting and ending mesh points on each
    // // processor for the interior part of the mesh.
    VecGetOwnershipRange(temp, &rstart, &rend);
    VecGetLocalSize(temp,&nlocal);

    Mat            A;                   // linear system matrix

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
    double start_time = MPI_Wtime();
    if (rstart == 0) 
    {
        rstart = 1;
        i      = 0; col[0] = 0; col[1] = 1; value[0] = 1 + 2.0*r1; value[1] = -r1;
        MatSetValues(A, 1, &i, 2, col, value, INSERT_VALUES);
    }

    if (rend == N) 
    {
        rend = N-1;
        i    = N-1; col[0] = N-2; col[1] = N-1; value[0] = -2.0*r1; value[1] = 1 + 2.0*r1;
        MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);
    }

    // Set entries corresponding to the mesh interior
    value[0] = -r1; value[1] = 1 + 2.0*r1; value[2] = -r1;
    for (i=rstart; i<rend; i++) 
    {
        col[0] = i-1; col[1] = i; col[2] = i+1;
        MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);
    }

    // Assemble the matrix
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    double end_time = MPI_Wtime();
    PetscPrintf(PETSC_COMM_WORLD, "Matrix assembly time: %.16f seconds\n", end_time - start_time);

    // MatView(A,PETSC_VIEWER_STDOUT_WORLD);

    initialize(temp, F, h5_tls, vec2arry); // initialize the temperature vector and source 

    // update_sourcevec(F, 1); // evaluate F at t = dt

    start_time = MPI_Wtime();
    time_loop_impliciteuler(temp, F, A, h5_tls, vec2arry);   // loop over the time step, meanwhile write the h5 file for each time step
    end_time = MPI_Wtime();
    PetscPrintf(PETSC_COMM_WORLD, "Sovler time: %.16f seconds\n", end_time - start_time);

    MatDestroy(&A);
}
