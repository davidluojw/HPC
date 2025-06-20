# HPC

## Compilation process (using CMake)

Please ensure that the following dependencies are installed:

- CMake
- Supported C/C++compilers (such as GCC)
- MPI libraries (MPICH)

### Compilation steps

1. Clone this repository:

```bash
git clone  https://github.com/davidluojw/HPC.git
cd HPC
```

2. Create a build directory and configure it:

```bash
mkdir build
cd build
cmake ..
```

>To specify an MPI compiler, set environment variables before the cmake command, for example:
>
> ```bash
> export CC=mpicc
> export CXX=mpicxx
> cmake ..
> ```

3. Compilation project:

```bash
make -j
```

After successful compilation, executable files are usually generated in the 'build/' directory.
