cmake_minimum_required(VERSION 3.12.0)
project(20250210 LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_executable(exercise-01 src/exercise-01.cpp src/ADR1D.cpp)
deal_ii_setup_target(exercise-01)

add_executable(exercise-02 src/exercise-02.cpp src/UStokes.cpp)
deal_ii_setup_target(exercise-02)