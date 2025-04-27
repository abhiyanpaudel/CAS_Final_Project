cmake_minimum_required(VERSION 3.10)
project(build_simplex_mesh)

find_package(Omega_h REQUIRED)

add_executable(build_simplex_mesh build_simplex_mesh.cpp)
target_link_libraries(build_simplex_mesh Omega_h::omega_h)
