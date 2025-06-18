#include "Vec2Array.hpp"
#include "petsc.h"

Vec2Array::Vec2Array(){}

std::vector<double> Vec2Array::get_vector_array(Vec vec){
    PetscScalar *array;
    PetscInt size;

    // get size
    VecGetSize(vec, &size);

    // gei array
    VecGetArray(vec, &array);

    // assign to array
    std::vector<double> data(size);

    for (PetscInt i = 0; i < size; i++){
        data[i] = static_cast<double>(array[i]);
    }

    // restore array
    VecRestoreArray(vec, &array);

    return data;
}


