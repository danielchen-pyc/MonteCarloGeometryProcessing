# MonteCarloGeometryProcessing

## Introduction
Course project for CSC2521 Topics in Computer Graphics: Physics-based Animation.

An implementation based on the paper by Rohan Sawhney and Keenan Crane \
["Monte Carlo Geometry Processing: A Grid-Free Approach to PDE-Based Methods on Volumetric Domains"](https://www.cs.cmu.edu/~kmcrane/Projects/MonteCarloGeometryProcessing/paper.pdf)

**This project was forked and build on top of [this repo](https://github.com/hgeorge21/MonteCarloGeometryProcessing). Most of the Monte Carlo Method are implemented from this code.**

## System Info
Tested on 
- 2018 MacBook Pro Sequoia (Intel i7)
- Apple clang version 16.0.0 (clang-1600.0.26.4)
- Target: x86_64-apple-darwin24.1.0
- `export CPATH="/Library/Developer/CommandLineTools/SDKs/MacOSX15.1.sdk/usr/include/c++/v1"`
- `export C_INCLUDE_PATH=/Library/Developer/CommandLineTools/SDKs/MacOSX15.1.sdk/usr/include/c++/v1`
- `export CPLUS_INCLUDE_PATH=/Library/Developer/CommandLineTools/SDKs/MacOSX15.1.sdk/usr/include/c++/v1`

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

## Building Note
- Follow the steps below (#Compilation)
- If cmake cannot find package igl::opengl or other packages:
    - Make sure `CMakeCache.txt` contains the correct ON/OFF values for each package
    - Make sure `add_subdirectory(${PROJECT_SOURCE_DIR}/libigl)` is presence in `CMakeLists.txt`
- If cannot find MacOS SDK files such as `cmath`
    - Set environment variable `export CPATH=/Library/Developer/CommandLineTools/SDKs/MacOSX15.0.sdk/usr/include/c++/v1`
- If `cmake` fails when installing `gmp`, remember to set `CGAL_WITH_GMPXX:BOOL=OFF`


## Program Structure
- `bary_rkpm.cpp` is the main program that runs the Monte-Carlo WoS algorithm to compute the Barycentric coordinates. `make` compiles this file into the `build/MonteCarloPDE` executable. 
- `test.cpp` creates the bounding box for a specified object in `data/` folder. Run `g++ -std=c++17 -I/usr/local/include/eigen3 -I./libigl/include/ -I./include/ test.cpp -o test` to compile the program into `test`. 
- `deliverables/` folder stores all the precompiled executables. See more info in #Precompiled Programs.
- `src/` folder stores all the helper functions and `include/` folder stores all the header files for these scripts.


## Precompiled Programs
Precompiled executables are in the `deliverables/` folder. Commands for running the solver will be printed in the terminal. 
- `sample_yk` program displays the sampled yk points (on the cage surface) after running the Monte Carlo solver.
- `bunny_bounding_box` program displays the bounding box for the bunny object
