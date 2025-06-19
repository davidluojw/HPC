#ifndef HEAT_SOLVER
#define HEAT_SOLVER

#include <petsc.h>
#include "Manufactured_Solution.hpp"
#include "hdf5_tools.hpp"
#include "Vec2Array.hpp"

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
    void initialize(Vec &temp, Vec &F, hdf5_tools * const & h5_tls, Vec2Array * const & vec2arry);

    // update temperature
    void update_temperature_expliciteuler(Vec &temp, Vec &F, Mat &A);

    // update source vec
    void update_sourcevec(Vec &F, int j);

    // loop time step
    void time_loop(Vec &temp, Vec &F, Mat &A, hdf5_tools * const & h5_tls, Vec2Array * const & vec2arry);


};

#endif