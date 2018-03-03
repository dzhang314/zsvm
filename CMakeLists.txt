cmake_minimum_required(VERSION 3.0)
project(zsvm)

set(CMAKE_CXX_STANDARD 17)

add_compile_options(
        -Wall -Wextra -pedantic -Werror -pedantic-errors
        -O3 -march=native -mtune=native)

add_executable(zsvm
        Particle.hpp
        Particle.cpp
        Permutation.hpp
        Permutation.cpp
        JacobiCoordinates.hpp
        JacobiCoordinates.cpp
        SphericalECGContext.hpp
        SphericalECGContext.cpp
        RealVariationalSolver.hpp
        RealVariationalSolver.cpp
        main.cpp)
