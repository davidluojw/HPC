// heat_solver.cpp
#include "heat_solver.hpp"
#include "hdf5_tools.hpp"
#include "Vec2Array.hpp"
#include <petsc.h>

heat_solver::heat_solver(Manufactured_Solution *mn_sol_in)
    : mn_sol(mn_sol_in) {}

heat_solver::~heat_solver() {}


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
        PetscPrintf(PETSC_COMM_SELF,"from rank %d , rstart = %d, rend = %d, i = %d, u_val = %f, x_coor = %f \n",rank, rstart, rend, i, u_val, mn_sol->x_coor[i+1]);
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

}

void heat_solver::update_temperature_expliciteuler(Vec &temp, Vec &F, Mat &A){
    Vec temp_np1;
    VecDuplicate(temp, &temp_np1);

    // temp_np1 = A*temp
    MatMult(A, temp, temp_np1);

    // temp_np1 = temp_np1 + F
    VecAXPY(temp_np1, 1.0, F);

    // Update temp
    VecCopy(temp_np1, temp);

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
        // std::cout << "x_coot: " << mn_sol->x_coor[i+1] << "tt: " << tt << ", r1 " << mn_sol->r1 << ", r2" << mn_sol->r2 << std::endl;
        // std:: cout << "f_val: " << f_val << std::endl;
        if (i == 0){
            f_val += mn_sol->r1 * mn_sol->g(mn_sol->x_coor[0], tt);
            // PetscPrintf(PETSC_COMM_SELF,"from rank %d , i = %d, f_val = %f\n",rank, i, f_val);
        }
        else if (i == mn_sol->N - 1){
            double h = mn_sol->r1 * 2 * mn_sol->dx / mn_sol->kappa * mn_sol->h(mn_sol->x_coor[mn_sol->N], tt, 1);
            // PetscPrintf(PETSC_COMM_SELF, "h: %.16f\n", h);
            f_val += h;
            // std:: cout << "h(x,t), f_val: " << f_val << std::endl;
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

    // update_temperature_expliciteuler(temp, F, A);

    for (int j = 1; j <= mn_sol->M; j++) {
        
        update_sourcevec(F, j);  // update the source vector F for the next time step
        update_temperature_expliciteuler(temp, F, A);

        // VecView(temp, PETSC_VIEWER_STDOUT_WORLD);
        // VecView(F, PETSC_VIEWER_STDOUT_WORLD);

        double tt = j * mn_sol->dt;

        std::vector<double> step_temp = vec2arry->get_vector_array(temp);
        step_temp.insert(step_temp.begin(), mn_sol->u(mn_sol->x_coor[0], tt));
        if (rank == 0){
            std::cout << "step_temp: "<< j << "\n";
            // for (int ii = 0; ii < step_temp.size(); ++ii){
            //     std::cout << step_temp[ii] << "\t";
            // }
            // std::cout << std::endl;

            h5_tls->write_hdf5(j, 0, step_temp);

        }
    }
    MPI_Barrier(PETSC_COMM_WORLD);


}


void heat_solver::explicitEuler(Vec &temp, Vec &F, hdf5_tools * const & h5_tls, Vec2Array * const & vec2arry){
    int N = mn_sol->N;
    double r1 = mn_sol->r1;
    PetscInt       i,col[3],rstart,rend,nlocal;
    PetscScalar    value[3];

    PetscMPIInt    rank,size;
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // // Identify the starting and ending mesh points on each
    // // processor for the interior part of the mesh.
    VecGetOwnershipRange(temp, &rstart, &rend);
    VecGetLocalSize(temp,&nlocal);

    Mat            A;                   // linear system matrix

    // // Create matrix.
    // // We pass in nlocal as the "local" size of the matrix to force it
    // // to have the same parallel layout as the vector created above.
    double start_time = MPI_Wtime();
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
        // rstart = 1;
        i      = 0; 
        col[0] = 0; col[1] = 1; 
        value[0] = 1 - 2.0*r1; value[1] = r1;
        MatSetValues(A, 1, &i, 2, col, value, INSERT_VALUES);
    }

    if (rend == N) 
    {
        // rend = N-1;
        i    = N-1; 
        col[0] = N-2; col[1] = N-1; 
        value[0] = 2.0*r1; value[1] = 1 - 2.0*r1;
        MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);
    }

    // Set entries corresponding to the mesh interior
    // value[0] = r1; value[1] = 1 - 2.0*r1; value[2] = r1;
    for (i = std::max(rstart, 1); i < std::min(rend, N-1); ++i) 
    {
        col[0] = i-1; col[1] = i; col[2] = i+1;
        value[0] = r1; value[1] = 1 - 2.0*r1; value[2] = r1;
        MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);
    }

    // Assemble the matrix
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    double end_time = MPI_Wtime();
    PetscPrintf(PETSC_COMM_WORLD, "Matrix Assembly time: %.16f seconds\n", end_time - start_time);

    MatView(A,PETSC_VIEWER_STDOUT_WORLD);

    initialize(temp, F, h5_tls, vec2arry); // initialize the temperature vector and source 

    double t0 = mn_sol->initial_time;

    std::vector<double> init_temp = vec2arry->get_vector_array(temp);  // get the array of Vec temp
    init_temp.insert(init_temp.begin(), mn_sol->u(mn_sol->x_coor[0], t0));
    if (rank == 0){
        std::cout << "init_temp: \n";
        for (int ii = 0; ii < init_temp.size(); ++ii){
            std::cout << init_temp[ii] << "\t";
        }
        std::cout << std::endl;

        h5_tls->write_hdf5(0, 0, init_temp);
    }
    MPI_Barrier(PETSC_COMM_WORLD);

    start_time = MPI_Wtime();
    time_loop_expliciteuler(temp, F, A, h5_tls, vec2arry);   // loop over the time step, meanwhile write the h5 file for each time step
    end_time = MPI_Wtime();
    PetscPrintf(PETSC_COMM_WORLD, "Sovler time: %.16f seconds\n", end_time - start_time);

     MatDestroy(&A);
}


void heat_solver::update_temperature_impliciteuler(Vec &temp, Vec &F, Mat &A){

    KSP ksp;
    PC pc;

    // F = temp + F
    // VecAXPY(F, 1.0, temp);
    VecAXPY(temp, 1.0, F);
    // PetscPrintf(PETSC_COMM_WORLD, "Updated temperature vector:\n");
    // VecView(temp, PETSC_VIEWER_STDOUT_WORLD);

     // Create the linear solver and set various options
    KSPCreate(PETSC_COMM_WORLD, &ksp);

    // Set operators. Here the matrix that defines the linear system
    // also serves as the preconditioning matrix.
    KSPSetOperators(ksp, A, A);
    KSPSetType(ksp, KSPGMRES);

    // Set linear solver defaults for this problem (optional).
    // - By extracting the KSP and PC contexts from the KSP context,
    //   we can then directly call any KSP and PC routines to set
    //   various options.
    // - The following four statements are optional; all of these
    //   parameters could alternatively be specified at runtime via
    //   KSPSetFromOptions();
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCJACOBI);
    KSPSetTolerances(ksp, 1.e-14, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

    // Set runtime options, e.g.,
    //     -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    // These options will override those specified above as long as
    // KSPSetFromOptions() is called _after_ any other customization routines.
    KSPSetFromOptions(ksp);

    // Solve the linear system
    KSPSolve(ksp,temp,temp);

    KSPDestroy(&ksp);

}

void heat_solver::time_loop_impliciteuler(Vec &temp, Vec &F, Mat &A, hdf5_tools * const & h5_tls, Vec2Array * const & vec2arry) {

    PetscInt rstart, rend;
    PetscMPIInt    rank,size, nlocal;
    VecGetOwnershipRange(F, &rstart, &rend);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    VecGetLocalSize(temp,&nlocal);

    for (int j = 0; j < mn_sol->M; j++) {
        VecSet(F, 0.0);  // reset F to zero before each time step
        update_sourcevec(F, j+1);    // F1, F2, F3, F4, .... 
        if (j == 1) 
        {
            PetscPrintf(PETSC_COMM_WORLD, "The first iteration F1: \n");
            VecView(F, PETSC_VIEWER_STDOUT_WORLD);
        }
        // VecView(F, PETSC_VIEWER_STDOUT_WORLD);
        update_temperature_impliciteuler(temp, F, A);  // U1, 2, U3, U4, ...
         if (j == 1)
        {
            PetscPrintf(PETSC_COMM_WORLD, "Solved temperature vector:\n");
            VecView(temp, PETSC_VIEWER_STDOUT_WORLD);
        }
        
        // VecView(temp, PETSC_VIEWER_STDOUT_WORLD);
        // VecView(F, PETSC_VIEWER_STDOUT_WORLD);

        double tt = (j+1) * mn_sol->dt;

        std::vector<double> step_temp = vec2arry->get_vector_array(temp);
        step_temp.insert(step_temp.begin(), mn_sol->u(mn_sol->x_coor[0], tt));
        if (rank == 0){
            // std::cout << "step_temp: "<< j+1<< "\n";
            // for (int ii = 0; ii < step_temp.size(); ++ii){
            //     std::cout << step_temp[ii] << "\t";
            // }
            // std::cout << std::endl;

            h5_tls->write_hdf5(j+1, 0, step_temp);

        }

        if (j == 1)
        {
            PetscPrintf(PETSC_COMM_WORLD, "Initial step completed at time %.2f\n", tt);
            PetscPrintf(PETSC_COMM_WORLD, "Solution vector after initial step:\n");
            VecView(temp, PETSC_VIEWER_STDOUT_WORLD);
        }

    }
    MPI_Barrier(PETSC_COMM_WORLD);

}

void heat_solver::implicitEuler(Vec &temp, Vec &F, hdf5_tools * const & h5_tls, Vec2Array * const & vec2arry){
    int N = mn_sol->N;
    double r1 = mn_sol->r1;
    PetscInt       i,col[3],rstart,rend,nlocal;
    PetscScalar    value[3];

    PetscMPIInt    rank,size;
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

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
        // rstart = 1;
        i      = 0; 
        col[0] = 0; col[1] = 1; 
        value[0] = 1 + 2.0*r1; value[1] = -r1;
        MatSetValues(A, 1, &i, 2, col, value, INSERT_VALUES);
    }

    if (rend == N) 
    {
        // rend = N-1;
        i    = N-1; 
        col[0] = N-2; col[1] = N-1; 
        value[0] = -2.0*r1; value[1] = 1 + 2.0*r1;
        MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);
    }

    // Set entries corresponding to the mesh interior
    // value[0] = -r1; value[1] = 1 + 2.0*r1; value[2] = -r1;
    for (i = std::max(rstart, 1); i < std::min(rend, N-1); ++i)
    {
        col[0] = i-1; col[1] = i; col[2] = i+1;
        value[0] = -r1; value[1] = 1 + 2.0*r1; value[2] = -r1;
        MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);
    }

    // Assemble the matrix
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    double end_time = MPI_Wtime();
    PetscPrintf(PETSC_COMM_WORLD, "Matrix assembly time: %.16f seconds\n", end_time - start_time);

    MatView(A,PETSC_VIEWER_STDOUT_WORLD);

    initialize(temp, F, h5_tls, vec2arry); // initialize the temperature vector and source, U0

    double t0 = mn_sol->initial_time;

    std::vector<double> init_temp = vec2arry->get_vector_array(temp);  // get the array of Vec temp
    init_temp.insert(init_temp.begin(), mn_sol->u(mn_sol->x_coor[0], t0));
    if (rank == 0){
        std::cout << "init_temp: \n";
        for (int ii = 0; ii < init_temp.size(); ++ii){
            std::cout << init_temp[ii] << "\t";
        }
        std::cout << std::endl;

        h5_tls->write_hdf5(0, 0, init_temp);
    }
    MPI_Barrier(PETSC_COMM_WORLD);
    VecView(temp,PETSC_VIEWER_STDOUT_WORLD);

    start_time = MPI_Wtime();
    time_loop_impliciteuler(temp, F, A, h5_tls, vec2arry);   // loop over the time step, meanwhile write the h5 file for each time step
    end_time = MPI_Wtime();
    PetscPrintf(PETSC_COMM_WORLD, "Sovler time: %.16f seconds\n", end_time - start_time);

    MatDestroy(&A);
}
