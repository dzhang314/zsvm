cmake_minimum_required(VERSION 3.9)
project(zsvm)

set(CMAKE_CXX_STANDARD 17)

add_compile_options(
        -Wall -Wextra -pedantic -Werror -pedantic-errors
        -O3 -march=native -mtune=native)

add_executable(zsvm
        Particle.hpp
        Permutation.hpp
        JacobiCoordinates.hpp
        RealSolver.hpp
        main.cpp)
