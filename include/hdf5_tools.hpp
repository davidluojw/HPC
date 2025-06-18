#ifndef HDF5_TOOLS
#define HDF5_TOOLS

#include "hdf5.h"
#include <iostream>
#include <cmath>
#include <vector>
#include "petsc.h"

class hdf5_tools
{
private:
    std::string filename;
    hid_t file_id;
public:
    // constructor
    hdf5_tools(const std::string& filename);

    // destructor
    ~hdf5_tools(){}; 

    // set up hdf5
    void setup_hdf5();

    // write hdf5
    void write_hdf5( const int &time_index, std::vector<double>& data, double time_value );


};


#endif