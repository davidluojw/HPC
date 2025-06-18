#ifndef HEAT_SOLVER
#define HEAT_SOLVER

#include <petsc.h>
#include "Manufactured_Solution.hpp"


class heat_solver
{
private:
    // physical parameters
    Manufactured_Solution *mn_sol;
public:
    // constructor
    heat_solver(Manufactured_Solution * mn_sol);
    
    // destructor
    ~heat_solver();

    // initialization
    void initialize(Vec &temp, Vec &F);


};

#endif