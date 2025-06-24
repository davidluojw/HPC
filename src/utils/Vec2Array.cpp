#include "Vec2Array.hpp"
#include "petsc.h"

Vec2Array::Vec2Array(){}

Vec2Array::~Vec2Array() {}

std::vector<double> Vec2Array::get_vector_array(Vec temp){
    PetscInt rstart, rend, local_size, global_size;
    PetscMPIInt    rank,size;
    PetscScalar *local_array;

    // geu MPI processor
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);


    // get global& local size
    VecGetSize(temp, &global_size);
    VecGetOwnershipRange(temp, &rstart, &rend);
    local_size = rend - rstart;

    // get local array
    VecGetArray(temp, &local_array);

    // assign to array
    std::vector<double> local_data(local_size);
    for(PetscInt i = 0; i < local_size; i++){
        local_data[i] = static_cast<double>(local_array[i]);
    }

    // restore array
    VecRestoreArray(temp, &local_array);

    // Collect chunks of data from all ranks in the communicator to the root rank (inverse of scatter)
    // Calculate the local size of each process and broadcast it to the main process
    std::vector<int> counts(size), displs(size);
    int local_count = static_cast<int>(local_size);
    // local_size: proc[0], 5; proc, 5 => counts[2]: 5 5
    // counts.data(): The function is to return a bare pointer (of type int * here) pointing to 
    //                the first element of the continuous storage area inside the vector,  = &counts[0]
    MPI_Gather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, PETSC_COMM_WORLD);

    // calculate the displacement
    std::vector<double> global_data;
    if (rank == 0) {
        displs[0] = 0;
        for (int i = 1; i < size; ++i) {
            displs[i] = displs[i - 1] + counts[i - 1];
        }
        // displs[2]: 0 5 
        global_data.resize(global_size);
    }

    // Collect all process data to the main process
    // local_data.data(): local data pointer
    // local_count: local data size
    // MPI_DOUBLE: data type
    // global_data.data(): global data pointer (valid only for main process)
    // counts.data(): the number of data per process
    // displs.data(): the displacement of each process data in the global array
    // 0: main process number
    MPI_Gatherv(local_data.data(), local_count, MPI_DOUBLE,
                global_data.data(), counts.data(), displs.data(), MPI_DOUBLE,
                0, PETSC_COMM_WORLD);
    

    return global_data;  // only process 0 return 
}


