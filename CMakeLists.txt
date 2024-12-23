cmake_minimum_required(VERSION 3.1)
project(MonteCarloPDE)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake) 

add_subdirectory(${PROJECT_SOURCE_DIR}/libigl)

### libIGL options: choose between header only and compiled static library
set(LIBIGL ${PROJECT_SOURCE_DIR}/libigl)
set(LIBIGL_INCLUDE_DIR ${LIBIGL}/include)
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" ON)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"                   ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"                     ON)

find_package(LIBIGL REQUIRED QUIET)

# Can install a fuller version of eigen @ https://gitlab.com/libeigen/eigen.git
set(EIGENROOT "${PROJECT_SOURCE_DIR}/libigl/external/eigen")
include_directories("${EIGENROOT}")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include/")
file(GLOB SRCFILES src/*.cpp)

# add_executable(${PROJECT_NAME} "main.cpp" ${SRCFILES})
add_executable(${PROJECT_NAME} "bary_rkpm.cpp" ${SRCFILES})
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl igl::glfw)