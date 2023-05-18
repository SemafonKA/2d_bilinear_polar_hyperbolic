cmake_minimum_required(VERSION 3.0.0)
project(
   2d_bilinear_polar_hyperbolic
   VERSION 0.1.0
   LANGUAGES CXX
)

add_library(
   IterSolvers
   STATIC
   ./lib/three_steppers/Resourses/IterSolvers.cpp
)

add_library(
   LU
   STATIC
   ./lib/three_steppers/Resourses/LU.cpp
)

add_library(
   SparseMatrix
   STATIC
   ./lib/three_steppers/Resourses/SparseMatrix.cpp
)

add_executable(${PROJECT_NAME} ./lib/main.cpp)

target_link_libraries(
   ${PROJECT_NAME}
   IterSolvers
   LU
   SparseMatrix
)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})


if(MSVC)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++20")
else()
   set(CMAKE_CXX_STANDARD 20)
endif(MSVC)

include(CPack)
