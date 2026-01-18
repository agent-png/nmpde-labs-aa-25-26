cmake_minimum_required(VERSION 3.12.0)
project(UStokes LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_executable(exercise-01 src/exercise-01.cpp src/UStokes.cpp)
deal_ii_setup_target(exercise-01)