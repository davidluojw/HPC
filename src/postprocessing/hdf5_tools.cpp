// hdf5_tools.cpp
#include "hdf5_tools.hpp"
#include <hdf5.h>

hdf5_tools::hdf5_tools(const std::string& filename_in): filename(filename_in){}


void hdf5_tools::setup_hdf5(){
    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    H5Fclose(file_id);
}

void hdf5_tools::write_hdf5(){


}