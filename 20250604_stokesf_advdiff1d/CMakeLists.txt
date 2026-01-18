cmake_minimum_required(VERSION 3.12.0)
project(tde2 LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_library(ARD1D STATIC src/ARD1D.cpp)
deal_ii_setup_target(ARD1D)

add_executable(exercise-01 src/exercise-01.cpp)
target_link_libraries(exercise-01 ARD1D)
deal_ii_setup_target(exercise-01)

add_executable(exercise-02 src/exercise-02.cpp src/Stokes.cpp)
deal_ii_setup_target(exercise-02)
