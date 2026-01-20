cmake_minimum_required(VERSION 3.12.0)
project(lab02 LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_library(Parabolic STATIC src/Parabolic.cpp)
deal_ii_setup_target(Parabolic)

add_executable(exercise-01 src/exercise-01.cpp)
target_link_libraries(exercise-01 Parabolic)
deal_ii_setup_target(exercise-01)

