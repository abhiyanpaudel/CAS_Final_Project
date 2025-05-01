cmake_minimum_required(VERSION 3.20)

project(point2meshMap LANGUAGES CXX)

enable_testing()
# Find required packages
find_package(meshfields REQUIRED)
find_package(Omega_h REQUIRED)
find_package(MPI REQUIRED)
find_package(pcms REQUIRED)
find_package(Catch2 REQUIRED)

# Direct PETSc configuration
set(PETSC_DIR "/lore/paudea/sources/petsc")
set(PETSC_ARCH "linux-gnu-gpu-kokkos")
set(PETSC_INCLUDES "${PETSC_DIR}/include" "${PETSC_DIR}/${PETSC_ARCH}/include")
set(PETSC_LIBRARIES "${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.so")

# Check if PETSc libraries exist
if(NOT EXISTS "${PETSC_LIBRARIES}")
  message(FATAL_ERROR "PETSc library not found at: ${PETSC_LIBRARIES}")
endif()

add_executable(point2MeshMap tests/test_point_to_mesh_map.cpp)
target_link_libraries(point2MeshMap PUBLIC
                    Omega_h::omega_h
                    pcms::pcms
					Catch2::Catch2
                    meshfields::meshfields
                    ${PETSC_LIBRARIES})

target_include_directories(point2MeshMap PUBLIC 
                          ${CMAKE_CURRENT_SOURCE_DIR}
                          ${PETSC_INCLUDES})

include(CTest)
include(Catch)
catch_discover_tests(point2MeshMap)
