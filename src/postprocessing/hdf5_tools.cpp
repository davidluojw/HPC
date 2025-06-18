// hdf5_tools.cpp
#include "hdf5_tools.hpp"
#include <hdf5.h>


hdf5_tools::hdf5_tools(const std::string& filename_in): filename(filename_in){}


void hdf5_tools::setup_hdf5(){
    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    // only write h5 file on main process
    if (rank == 0){
        hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if(file_id < 0) {
            std::cerr << "Error creating HDF5 file: " << filename << std::endl;
            return;
        }

        herr_t status = H5Fclose(file_id);
        if(status < 0) {
            std::cerr << "Error closing HDF5 file: " << filename << std::endl;
        }
    }

    
}

void hdf5_tools::write_hdf5( const int &time_index, std::vector<double>& data, double time_value ){

    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // only write h5 file on main process
    if (rank == 0){
        hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        if(file_id < 0) {
            // 文件不存在，创建它
            file_id = H5Fcreate(filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
            if(file_id < 0) {
                std::cerr << "Error creating HDF5 file: " << filename << std::endl;
                return;
            }
        }

        const std::string group_name = std::to_string(900000000 + time_index);

        hid_t group_id = H5Gcreate(file_id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if(group_id < 0) {
            std::cerr << "Error creating group: " << group_name << std::endl;
            H5Fclose(file_id);
            return;
        }

        // Now write the data into the h5 file
        hsize_t dims[1] = { static_cast<hsize_t>(data.size()) };
        if(dims[0] > 0)
        {
            hid_t dataspace = H5Screate_simple(1, dims, NULL);
            if(dataspace < 0) {
                std::cerr << "Error creating dataspace" << std::endl;
                H5Gclose(group_id);
                H5Fclose(file_id);
                return;
            }
            hid_t dataset   = H5Dcreate( group_id, std::to_string(time_index).c_str(), H5T_NATIVE_DOUBLE, 
                                        dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
            if(dataset < 0) {
                std::cerr << "Error creating dataset" << std::endl;
                H5Sclose(dataspace);
                H5Gclose(group_id);
                H5Fclose(file_id);
                return;
            }

            herr_t status = H5Dwrite( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0] );
            if(status < 0) {
                std::cerr << "Error writing dataset" << std::endl;
            }

            // check_error(status, "write_doubleVector");
            H5Dclose( dataset );
            H5Sclose( dataspace );
        }   

        H5Gclose(group_id); 
        H5Fclose(file_id);
    }

    

}