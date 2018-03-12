cmake_minimum_required(VERSION 3.0)
project(zsvm)

set(CMAKE_CXX_STANDARD 17)

add_definitions(-DEIGEN_NO_DEBUG=1)

add_compile_options(
        -Wall -Wextra -pedantic -Werror -pedantic-errors
        -O3 -march=native -mtune=native)

set_source_files_properties(PackedLinearAlgebra.cpp
        PROPERTIES COMPILE_FLAGS -ffast-math)

add_executable(zsvm
        Particle.hpp
        Particle.cpp
        Permutation.hpp
        Permutation.cpp
        JacobiCoordinates.hpp
        JacobiCoordinates.cpp
        PackedLinearAlgebra.hpp
        PackedLinearAlgebra.cpp
        SphericalECGContext.hpp
        SphericalECGContext.cpp
        RealVariationalSolver.hpp
        RealVariationalSolver.cpp
        AmoebaOptimizer.hpp
        AmoebaOptimizer.cpp
        ScriptToken.hpp
        ScriptToken.cpp
        ScriptTokenizer.hpp
        ScriptTokenizer.cpp
        main.cpp)
