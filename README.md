# MonteCarloGeometryProcessing

## Note
- Follow the steps below (#Compilation)
- If cmake cannot find package igl::opengl or other packages:
    - Make sure `CMakeCache.txt` contains the correct ON/OFF values for each package
    - Make sure `add_subdirectory(${PROJECT_SOURCE_DIR}/libigl)` is presence in `CMakeLists.txt`
- If cannot find MacOS SDK files such as `cmath`
    - Set environment variable `export CPATH=/Library/Developer/CommandLineTools/SDKs/MacOSX15.0.sdk/usr/include/c++/v1`
- If `cmake` fails when installing `gmp`, remember to set `CGAL_WITH_GMPXX:BOOL=OFF`



## Introduction
Course project for [CSC2520 Geometry Processing](https://github.com/alecjacobson/geometry-processing-csc2520).

An implementation based on the paper by Rohan Sawhney and Keenan Crane \
["Monte Carlo Geometry Processing: A Grid-Free Approach to PDE-Based Methods on Volumetric Domains"](https://www.cs.cmu.edu/~kmcrane/Projects/MonteCarloGeometryProcessing/paper.pdf)

## Compilation
This project makes use of open-source libraries `libigl` and `Eigen` which will be included by cloning `libigl` repository.

To clone the repository, use 
```
git clone https://github.com/hgeorge21/MonteCarloGeometryProcessing.git
cd MonteCarloGeometryProcessing
git clone --recursive https://github.com/libigl/libigl.git
```

To build the project, go into the project directory and do the following
```
mkdir build
cd build
cmake ..
make    # to compile the program to MonteCarloPDE
```
