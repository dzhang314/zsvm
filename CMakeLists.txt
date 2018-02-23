cmake_minimum_required(VERSION 3.9)
project(zsvm)

set(CMAKE_CXX_STANDARD 17)

add_executable(zsvm Particle.hpp JacobiCoordinates.hpp main.cpp RealSolver.hpp)
