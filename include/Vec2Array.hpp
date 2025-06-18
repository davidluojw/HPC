#ifndef VEC2ARRAY
#define VEC2ARRAY

#include "petsc.h"
#include <iostream>
#include <cmath>
#include <vector>

class Vec2Array
{
private:
    /* data */
public:
    Vec2Array();
    ~Vec2Array();

    std::vector<double> get_vector_array(Vec temp);
};




#endif