cmake_minimum_required(VERSION 3.8)
project(zsvm)

set(CMAKE_CXX_STANDARD 17)

add_definitions(-DEIGEN_NO_DEBUG=1)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	include_directories(SYSTEM "/mnt/c/Programs/boost_1_67_0")
	include_directories(SYSTEM "/mnt/c/Programs/eigen-eigen-5a0156e40feb")
    add_compile_options(
        -Wall -Wextra -pedantic -Werror -pedantic-errors
        -O3 -march=native -mtune=native -fopenmp)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # using Intel C++
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	include_directories(SYSTEM "C:/Programs/boost_1_67_0")
	include_directories(SYSTEM "C:/Programs/eigen-eigen-5a0156e40feb")
endif()

add_executable(zsvm
        EnumTypes.hpp
        Particle.hpp
        Permutation.hpp
        JacobiCoordinates.hpp
        PackedLinearAlgebra.hpp
        PackedLinearAlgebraImpl.hpp
        RealVariationalSolver.hpp
        SphericalECGDispersionContext.hpp
        SphericalECGConfiningContext.hpp
        SphericalECGPairwiseContext.hpp
        SphericalECGOverlapContext.hpp
        SphericalECGJacobiContext.hpp
        SphericalECGJacobiVariationalOptimizer.hpp
        AmoebaOptimizer.hpp
        ScriptToken.hpp
        ScriptToken.cpp
        ScriptTokenizer.hpp
        ScriptTokenizer.cpp
        ScriptCommand.hpp
        ScriptCommand.cpp
        ScriptParser.hpp
        ScriptParser.cpp
        main.cpp)

find_package(OpenMP)
target_link_libraries(zsvm ${OpenMP_CXX_LIBRARIES} ${PROJECT_LINK_LIBS})
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
