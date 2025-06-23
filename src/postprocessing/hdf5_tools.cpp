// hdf5_tools.cpp
#include "hdf5_tools.hpp"
#include <hdf5.h>
#include <iostream>
#include <cmath>
#include <vector>


hdf5_tools::hdf5_tools(const std::string& filename_in): filename(filename_in){}

hdf5_tools::~hdf5_tools() {}


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

void hdf5_tools::write_hdf5( const int &time_index, const int &dataset_name,  std::vector<double>& data){

    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // only read h5 file on main process
    if (rank == 0){
        hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        if(file_id < 0) {
            // File does not exist, create it
            // file_id = H5Fcreate(filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
            // if(file_id < 0) {
                // std::cerr << "Error creating HDF5 file: " << filename << std::endl;
                // return;
            // }
            std::cerr << "Error opening HDF5 file for writing: " << filename << std::endl;
            return;
        }

        const std::string group_name = std::to_string(900000000 + time_index);

        // std::cout << "Writing group: " << group_name << std::endl;

        // 1. Define the original error handling function pointer and client data
        H5E_auto2_t old_func;
        void *old_client_data;
        // 2. Turn off error output before attempting to open the group
        H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
        H5Eset_auto(H5E_DEFAULT, NULL, NULL);

        // 3. Attempt to open the group (error will not print at this time)
        hid_t group_id = H5Gopen(file_id, group_name.c_str(), H5P_DEFAULT);

        // 4. Immediately restore error handling settings
        H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);

        // 5. Check if the group exists and create it
        if (group_id < 0) {
            // Group does not exist, create it
            group_id = H5Gcreate(file_id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if(group_id < 0) {
                std::cerr << "Error creating group: " << group_name << std::endl;
                H5Fclose(file_id);
                return;
            }
        }

        // hid_t group_id = H5Gcreate(file_id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        // if(group_id < 0) {
        //     std::cerr << "Error creating group: " << group_name << std::endl;
        //     H5Fclose(file_id);
        //     return;
        // }

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
            hid_t dataset   = H5Dcreate( group_id, std::to_string(dataset_name).c_str(), H5T_NATIVE_DOUBLE, 
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


void hdf5_tools::read_h5(const std::string& h5file_name, std::vector<std::vector<double>>& data){
    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // only write h5 file on main process
    if (rank == 0){
        hid_t fid = H5Fopen(h5file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);  // H5F_ACC_RDONLY: read only
        if (fid < 0) std::cerr << "Cannot open HDF5 file" << std::endl;

        // 1. Collect group name (loop time step)
        hsize_t nobj;  
        H5Gget_num_objs(fid, &nobj);  //H5Gget_num_objs: obtain the number of objects in the root directory
        std::vector<std::string> groups;
        for (hsize_t i = 0; i < nobj; ++i)
        {
            char name[64]; 
            H5Gget_objname_by_idx(fid, i, name, 64); // obtain group names by index, i.e. "900000000"
            groups.emplace_back(name);  // insert a new element constructed with name as a parameter at the end of the groups container
        }

        // for (size_t i = 0; i < groups.size(); ++i) {
        //     std::cout << "Found group: " << groups[i] << std::endl;
        // }

        // 2. Data size
        hid_t g0 = H5Gopen(fid, groups[0].c_str(), H5P_DEFAULT);  // Open the first time step group
        hid_t d0 = H5Dopen(g0, "0", H5P_DEFAULT);                 // Open dataset '0', temperature field
        hid_t sp = H5Dget_space(d0);                              // Obtain dataset space
        hsize_t nx;                  
        H5Sget_simple_extent_dims(sp, &nx, nullptr);              // Read data length (nx)
        H5Dclose(d0); 
        H5Gclose(g0);

        const int nt = static_cast<int>(groups.size());           // groups size,  size_t transfer to int
        // Initialize 2D array: [time steps x data length],
        // Create a two-dimensional array of data based on the size of the groups
        // assign():Set data to have nt elements, each element is a std:: vector<double>of size nx, and all values are initialized to 0.0 by default
        // Therefore, data is a 2D array
        data.assign(nt, std::vector<double>(nx));

        // 3. Read all (t, x)
        for (int it = 0; it < nt; ++it)
        {
            hid_t g  = H5Gopen(fid, groups[it].c_str(), H5P_DEFAULT);  // Open time step group
            hid_t ds = H5Dopen(g, "0", H5P_DEFAULT);                   // Open the dataset "0", temperature field
            H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, data[it].data());                    // Read dataset ds to data array
            H5Dclose(ds); 
            H5Gclose(g);
        }
        H5Fclose(fid);
    }
}