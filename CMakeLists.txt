cmake_minimum_required(VERSION 3.10)

# Project name and language
project(final_project)

# 指定C++标准
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set the path to your PETSc installation
set(HOME_DIR /Users/lawkawai) 
set(VTK_DIR /Users/lawkawai/lib/VTK-9.3.1-shared/lib/cmake/vtk-9.3)
set(PETSC_DIR ${HOME_DIR}/lib/petsc-3.22.1-opt-h)
set(PETSC_ARCH .)
set(HDF5_ROOT ${PETSC_DIR})


set(CMAKE_C_COMPILER ${HOME_DIR}/lib/mpich-4.2.3-opt/bin/mpicc)
set(CMAKE_CXX_COMPILER ${HOME_DIR}/lib/mpich-4.2.3-opt/bin/mpicxx) 
set(CMAKE_CXX_FLAGS "-O0 -W -Wshadow -Wall -Wextra -Wuninitialized -Wno-unused-parameter")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MACOSX_RPATH 1)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

# Find the PETSc package
find_package(VTK REQUIRED)
find_package(PETSc REQUIRED)
find_package(HDF5 REQUIRED)

include_directories(${VTK_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${PETSC_INC})
include_directories(./include)

add_library(solver_heat ./src/solver/heat_solver.cpp)
add_library(post_hdf5 ./src/postprocessing/hdf5_tools.cpp)
add_library(vec2arry ./src/utils/Vec2Array.cpp)
target_link_libraries(solver_heat ${PETSC_LIB})
target_link_libraries(post_hdf5 ${PETSC_LIB})
target_link_libraries(vec2arry ${PETSC_LIB})
# set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${VTK_LIBRARIES})
# set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${HDF5_LIBRARIES})
# set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_LIB})
# message(STATUS "External Libraries: " ${EXTRA_LINK_LIBS})

# Add an executable for the project
add_executable(driver_heat driver_1D_transient_heat.cpp)

# Link PETSc to the project
target_link_libraries(driver_heat ${PETSC_LIB} ${HDF5_LIBRARIES} ${VTK_LIBRARIES})
target_link_libraries(driver_heat solver_heat)
target_link_libraries(driver_heat post_hdf5)
target_link_libraries(driver_heat vec2arry)
# Link Libraries


# Compiler flags and options
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PETSC_C_FLAGS}")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PETSC_CXX_FLAGS}")
