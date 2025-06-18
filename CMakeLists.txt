cmake_minimum_required(VERSION 3.10)

# Project name and language
project(final_project)

# 指定C++标准
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set the path to your PETSc installation
set(HOME_DIR /Users/lawkawai) 
set(PETSC_DIR ${HOME_DIR}/lib/petsc-3.22.1-opt-h)
set(PETSC_ARCH .)

set(CMAKE_C_COMPILER ${HOME_DIR}/lib/mpich-4.2.3-opt/bin/mpicc)
set(CMAKE_CXX_COMPILER ${HOME_DIR}/lib/mpich-4.2.3-opt/bin/mpicxx) 
set(CMAKE_CXX_FLAGS "-O0 -W -Wshadow -Wall -Wextra -Wuninitialized -Wno-unused-parameter")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MACOSX_RPATH 1)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

# Find the PETSc package
find_package(PETSc REQUIRED)

include_directories(${PETSC_INC})
include_directories(./include)

add_library(solver_heat ./src/solver/heat_solver.cpp)
target_link_libraries(solver_heat ${PETSC_LIB})
# set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_LIB})
# message(STATUS "External Libraries: " ${EXTRA_LINK_LIBS})

# Add an executable for the project
add_executable(driver_heat driver_1D_transient_heat.cpp)

# Link PETSc to the project
target_link_libraries(driver_heat ${PETSC_LIB})
target_link_libraries(driver_heat solver_heat)


# Link Libraries


# Compiler flags and options
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PETSC_C_FLAGS}")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PETSC_CXX_FLAGS}")
