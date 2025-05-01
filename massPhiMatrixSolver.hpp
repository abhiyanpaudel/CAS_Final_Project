/**
 * @file massPhiMatrixSolver.hpp
 * @brief Solver for the Mass-Phi linear system in particle-to-mesh mapping
 * @author [Author Name]
 * @date May 1, 2025
 *
 * This file contains functions to solve the linear system Mass * x = Phi * source_values
 * for mapping field values from source points to a target mesh using a finite
 * element approach with PETSc and Kokkos for GPU acceleration.
 */

#ifndef MASS_PHI_SOLVER_HPP
#define MASS_PHI_SOLVER_HPP

#include <pcms/point_search.h>
#include <Omega_h_shape.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_for.hpp>
#include <MeshField_Shape.hpp>
#include <petscmat_kokkos.hpp>
#include <petscvec_kokkos.hpp>
#include <petscvec.h>
#include <petscsys.h>
#include <petscksp.h>
#include <Omega_h_array.hpp>
#include <Kokkos_Core.hpp>

#include "calculateMassMatrix.hpp"
#include "calculatePhiMatrix.hpp"

/**
 * @brief Creates a PETSc vector with Kokkos backend from Omega_h data
 *
 * Converts Omega_h::Reals data into a PETSc vector that can be used with
 * Kokkos-enabled PETSc operations for GPU acceleration.
 *
 * @param len Length of the vector to create
 * @param data Source Omega_h::Reals data
 * @return Vec PETSc vector initialized with the data
 */
static Vec createKokkosVec(PetscInt len, const Omega_h::Reals &data) {
  std::cout << "DEBUG: createKokkosVec - Creating vector of length " << len << std::endl;
  Vec v;
  PetscErrorCode ierr;
  ierr = VecCreateSeqKokkos(PETSC_COMM_SELF, len, &v);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  
  std::cout << "DEBUG: createKokkosVec - Vector created, now filling with data" << std::endl;
  PetscScalar  *array;
  ierr = VecGetArray(v, &array);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  
  auto host_data = Omega_h::HostRead(data);
  
  for (PetscInt i = 0; i < len; ++i) {
	  array[i] = host_data[i];
  }
  
  ierr = VecRestoreArray(v, &array);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  std::cout << "DEBUG: createKokkosVec - Vector filled and ready" << std::endl;
  return v;
}

/**
 * @brief Performs matrix-vector multiplication using PETSc and Kokkos
 *
 * Multiplies the matrix M by the vector u to produce a new vector.
 * Operations are performed on GPU if available through Kokkos.
 *
 * @param M The matrix
 * @param u The vector to multiply with
 * @return Vec Result of M * u
 */
static Vec multiplyMatVec(Mat M, Vec u) {
  std::cout << "DEBUG: multiplyMatVec - Start" << std::endl;
  PetscInt m,n;
  PetscErrorCode ierr;

  ierr = MatGetSize(M, &m, &n);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  std::cout << "DEBUG: multiplyMatVec - Matrix size: " << m << " x " << n << std::endl;
  
  Vec b;
  ierr = VecCreateSeqKokkos(PETSC_COMM_SELF, m, &b);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  
  std::cout << "DEBUG: multiplyMatVec - Performing matrix-vector multiplication" << std::endl;
  ierr = MatMult(M, u, b);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  std::cout << "DEBUG: multiplyMatVec - Complete" << std::endl;

  return b;
}

/**
 * @brief Solves a linear system Ax = b using PETSc's KSP solvers
 *
 * Uses PETSc's Krylov Subspace solvers to find x in Ax = b.
 * The solver can be configured through PETSc runtime options.
 *
 * @param A The system matrix
 * @param b The right-hand side vector
 * @return Vec Solution vector x
 */
static Vec solveLinearSystem(Mat A, Vec b) {
  std::cout << "DEBUG: solveLinearSystem - Start" << std::endl;
  PetscInt m,n;
  PetscErrorCode ierr;
  
  ierr = MatGetSize(A, &m, &n);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  std::cout << "DEBUG: solveLinearSystem - Matrix size: " << m << " x " << n << std::endl;
  
  Vec x;
  ierr = VecCreateSeqKokkos(PETSC_COMM_SELF, n, &x);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  
  KSP ksp;
  ierr = KSPCreate(PETSC_COMM_SELF, &ksp);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  
  ierr = KSPSetOperators(ksp, A, A);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
 
  ierr = KSPSetFromOptions(ksp);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  
  ierr = KSPSetUp(ksp);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  
  std::cout << "DEBUG: solveLinearSystem - Solving system" << std::endl;
  ierr = KSPSolve(ksp, b, x);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  std::cout << "DEBUG: solveLinearSystem - System solved" << std::endl;

  ierr = KSPDestroy(&ksp);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  
  return x;
}

/**
 * @brief Constructs the mass matrix for the target mesh
 *
 * Creates a mass matrix using the calculateMassMatrix function.
 * The mass matrix represents the inner product of basis functions.
 *
 * @param mesh The target Omega_h mesh
 * @return Mat PETSc matrix representing the mass matrix
 */
static Mat buildMassMatrixKokkos(Omega_h::Mesh& mesh){
  std::cout << "DEBUG: buildMassMatrixKokkos - Building mass matrix" << std::endl;
  Mat mass;
  PetscErrorCode ierr = calculateMassMatrix(mesh, &mass);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  std::cout << "DEBUG: buildMassMatrixKokkos - Mass matrix built" << std::endl;
  return mass;
}

/**
 * @brief Constructs the Phi matrix for mapping from source points to target mesh
 *
 * Creates the Phi matrix using the calculatePhiMatrix function.
 * The Phi matrix defines how source point values are mapped to mesh vertices.
 *
 * @param mesh The target Omega_h mesh
 * @param src_coordinates Coordinates of the source points
 * @return Mat PETSc matrix representing the Phi mapping matrix
 */
static Mat buildPhiMatrixKokkos(Omega_h::Mesh& mesh,
								const Omega_h::Reals& src_coordinates){
  std::cout << "DEBUG: buildPhiMatrixKokkos - Building phi matrix" << std::endl;
  Mat phi;
  PetscErrorCode ierr = calculatePhiMatrix(mesh, src_coordinates, &phi);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  std::cout << "DEBUG: buildPhiMatrixKokkos - Phi matrix built" << std::endl;
  return phi;
}

/**
 * @brief Main solver function that maps source point values to a target mesh
 *
 * Solves the linear system Mass * x = Phi * source_values to find the 
 * coefficients that represent the best mapping from source point values
 * to the target mesh in a finite element sense.
 * 
 * Algorithm:
 * 1. Build the mass and Phi matrices
 * 2. Compute right-hand side b = Phi * source_values
 * 3. Solve Mass * x = b for x
 * 4. Return x as the solution vector
 *
 * @param mesh The target Omega_h mesh
 * @param source_coordinates Coordinates of the source points
 * @param source_values Values at the source points
 * @return Omega_h::Reals Solution vector representing mapped values at mesh vertices
 * @throws std::runtime_error If source_values size doesn't match source_coordinates
 */
Omega_h::Reals solveMassPhiSystem(
  Omega_h::Mesh &mesh,
  const Omega_h::Reals &source_coordinates,
  const Omega_h::Reals &source_values
) {
  std::cout << "DEBUG: solveMassPhiSystem - Start" << std::endl;
  PetscInt ncols = source_coordinates.size()/mesh.dim();
  std::cout << "DEBUG: solveMassPhiSystem - Source points: " << ncols << std::endl;
  
  if ((PetscInt)source_values.size() != ncols) {
    std::cerr << "ERROR: source_values size (" << source_values.size() 
              << ") doesn't match expected size (" << ncols << ")" << std::endl;
    throw std::runtime_error("source_values length mismatch");
  }

  std::cout << "DEBUG: solveMassPhiSystem - Building matrices" << std::endl;
  Mat mass = buildMassMatrixKokkos(mesh);
  Mat phi = buildPhiMatrixKokkos(mesh, source_coordinates);

  std::cout << "DEBUG: solveMassPhiSystem - Creating vectors and solving system" << std::endl;
  Vec u = createKokkosVec(ncols, source_values);
  Vec b = multiplyMatVec(phi, u);
  Vec x = solveLinearSystem(mass, b);

  std::cout << "DEBUG: solveMassPhiSystem - Extracting solution" << std::endl;
  Omega_h::Write<Omega_h::Real> solution_vector(ncols, 0, "stores the solution coefficients");
  {
    PetscErrorCode ierr;
    PetscScalar* array;
    ierr = VecGetArray(x, &array);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    auto solution_host = Omega_h::HostWrite<Omega_h::Real>(ncols);

    for (PetscInt i = 0; i < ncols; ++i){
      solution_host[i] = array[i];
    }

    solution_vector = Omega_h::Write<Omega_h::Real>(solution_host);
    ierr = VecRestoreArray(x, &array);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
  }

  std::cout << "DEBUG: solveMassPhiSystem - Cleaning up" << std::endl;
  PetscErrorCode ierr;
  ierr = VecDestroy(&u);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  
  ierr = VecDestroy(&b);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  
  ierr = VecDestroy(&x);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  
  ierr = MatDestroy(&phi);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  
  ierr = MatDestroy(&mass);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  
  std::cout << "DEBUG: solveMassPhiSystem - Complete" << std::endl;
  return Omega_h::read(solution_vector);
}

#endif // MASS_PHI_SOLVER_HPP
