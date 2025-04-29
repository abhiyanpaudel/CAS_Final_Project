#include "calculateMassMatrix.hpp"
#include "massPhiMatrixSolver.hpp"
#include <Omega_h_build.hpp>
#include <Omega_h_tag.hpp>
#include <petsc.h>
#include <petscvec_kokkos.hpp>
#include <petscmat_kokkos.hpp>
#include <fenv.h>

Omega_h::Reals evalFuncValues(const Omega_h::Reals& coordinates, const int dim){
	std::cout << "DEBUG: evalFuncValues - Start" << std::endl;
	int num_points = coordinates.size()/dim; 
	std::cout << "DEBUG: evalFuncValues - Number of points: " << num_points << std::endl;
	
	Omega_h::Write<Omega_h::Real> funcValues(num_points, 0.0, "stores the function values at the given coordinates");
	Omega_h::parallel_for(num_points, OMEGA_H_LAMBDA(const int i){
		auto x = coordinates[i * dim];
		auto y = coordinates[i * dim + 1];
		funcValues[i] = x  + 3 * y;  
	});
	Kokkos::fence();
	
	std::cout << "DEBUG: evalFuncValues - Completed" << std::endl;
	return Omega_h::read(funcValues);
}

int main(int argc, char** argv) {
  std::cout << "DEBUG: Program start" << std::endl;

  Kokkos::ScopeGuard kokkosGuard(argc, argv);
  auto lib = Omega_h::Library(&argc, &argv);

  std::cout << "DEBUG: Kokkos and Omega_h initialized" << std::endl;

  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  std::cout << "DEBUG: PETSc initialized" << std::endl;

  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " sourceMesh.osh targetMesh.osh\n";
    return EXIT_FAILURE;
  }

  std::cout << "DEBUG: Loading source mesh from: " << argv[1] << std::endl;
  Omega_h::Mesh srcMesh(&lib);
  Omega_h::binary::read(argv[1], lib.world(), &srcMesh);
  std::cout << "DEBUG: Source mesh loaded. Elements: " << srcMesh.nelems() << ", Vertices: " << srcMesh.nverts() << std::endl;
  
  std::cout << "DEBUG: Loading target mesh from: " << argv[2] << std::endl;
  Omega_h::Mesh tgtMesh(&lib);
  Omega_h::binary::read(argv[2], lib.world(), &tgtMesh);
  std::cout << "DEBUG: Target mesh loaded. Elements: " << tgtMesh.nelems() << ", Vertices: " << tgtMesh.nverts() << std::endl;

  auto source_coordinates = srcMesh.coords();
  std::cout << "DEBUG: Source coordinates array size: " << source_coordinates.size() << std::endl;
  
  auto source_values = evalFuncValues(source_coordinates, srcMesh.dim());
  std::cout << "DEBUG: Source values computed. Array size: " << source_values.size() << std::endl;

  std::cout << "DEBUG: Calling solveMassPhiSystem" << std::endl;
  Omega_h::Reals solution_vector = solveMassPhiSystem(
    tgtMesh,
    source_coordinates,
    source_values
  );
  std::cout << "DEBUG: solveMassPhiSystem completed. Solution vector size: " << solution_vector.size() << std::endl;

  std::cout << "DEBUG: Adding solution tag to target mesh" << std::endl;
  tgtMesh.add_tag(0,"solution",1, solution_vector);
  
  std::cout << "DEBUG: Writing solution to VTK file" << std::endl;
  Omega_h::vtk::write_parallel(
    "solution.vtk", &tgtMesh, 0
  );
  std::cout << "DEBUG: VTK file written" << std::endl;

  std::cout << "DEBUG: Finalizing PETSc" << std::endl;
  PetscCall(PetscFinalize());
  std::cout << "DEBUG: Program completed successfully" << std::endl;
  return EXIT_SUCCESS;
}
