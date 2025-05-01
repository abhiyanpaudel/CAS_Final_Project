#include "calculateMassMatrix.hpp"
#include "massPhiMatrixSolver.hpp"
#include <Omega_h_build.hpp>
#include <Omega_h_tag.hpp>
#include <petsc.h>
#include <petscvec_kokkos.hpp>
#include <petscmat_kokkos.hpp>
#include <fenv.h>
#include <string>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_session.hpp>

Omega_h::Reals evalFuncValues(const Omega_h::Reals& coordinates, const int dim, const int deg){
	std::cout << "DEBUG: evalFuncValues - Start" << std::endl;
	int num_points = coordinates.size()/dim; 
	std::cout << "DEBUG: evalFuncValues - Number of points: " << num_points << std::endl;
	
	Omega_h::Write<Omega_h::Real> funcValues(num_points, 0.0, "stores the function values at the given coordinates");
	Omega_h::parallel_for(num_points, OMEGA_H_LAMBDA(const int i){
		auto x = coordinates[i * dim];
		auto y = coordinates[i * dim + 1];
		if (deg == 1) {
		  funcValues[i] = x  + 3 * y;  
		} else { 
		   funcValues[i] = 3.0; // constant function
		}
	});
	Kokkos::fence();
	
	std::cout << "DEBUG: evalFuncValues - Completed" << std::endl;
	return Omega_h::read(funcValues);
}

static const std::string SRC_PATH = "/lore/paudea/projects/particle2mesh_map/create_mesh/source_mesh/source_mesh.osh";
static const std::string TGT_PATH = "/lore/paudea/projects/particle2mesh_map/create_mesh/target_mesh/target_mesh.osh";

static void read_mesh(Omega_h::Mesh& mesh, Omega_h::Library& lib, const std::string& path){

	Omega_h::binary::read(path.c_str(), lib.world(), &mesh);
}

TEST_CASE("constant function mapping", "[massphi][deg0]") {

  Omega_h::Library lib;

  Omega_h::Mesh src_mesh(&lib);
  Omega_h::Mesh tgt_mesh(&lib);

  read_mesh(src_mesh, lib, SRC_PATH);
  read_mesh(tgt_mesh, lib, TGT_PATH);

  auto src_coords = src_mesh.coords();
  auto tgt_coords = tgt_mesh.coords();
  REQUIRE(src_coords.size() > 0);
  REQUIRE(tgt_coords.size() > 0);

  // Degree 0 => constant
  auto src_vals = evalFuncValues(src_coords, 0);
  auto exact    = evalFuncValues(tgt_coords, 0);

  auto sol = solveMassPhiSystem(tgt_mesh, src_coords, src_vals);

  auto host_sol = Omega_h::HostRead(sol);
  auto host_exact = Omega_h::HostRead(exact);
  REQUIRE(sol.size() == exact.size());
  double tol = 1e-8;
  for (std::size_t i = 0; i < host_sol.size(); ++i) {
    CAPTURE(i, host_exact[i], host_sol[i]);
    CHECK(host_sol[i] == Approx(host_exact[i]).margin(tol));
  }

}

TEST_CASE("linear function mapping", "[massphi][deg1]") {

  Omega_h::Library lib;

  Omega_h::Mesh src_mesh(&lib);
  Omega_h::Mesh tgt_mesh(&lib);

  read_mesh(src_mesh, lib, SRC_PATH);
  read_mesh(tgt_mesh, lib, TGT_PATH);

  auto src_coords = src_mesh.coords();
  auto tgt_coords = tgt_mesh.coords();

  // Degree 1 => linear
  auto src_vals = evalFuncValues(src_coords, 1);
  auto exact    = evalFuncValues(tgt_coords, 1);

  auto sol = solveMassPhiSystem(tgt_mesh, src_coords, src_vals);

  auto host_sol = Omega_h::HostRead(sol);
  auto host_exact = Omega_h::HostRead(exact);
  REQUIRE(sol.size() == exact.size());
  double tol = 1e-7;
  for (std::size_t i = 0; i < sol.size(); ++i) {
    CAPTURE(i, exact[i], sol[i]);
    CHECK(host_sol[i] == Approx(host_exact[i]).margin(tol));
  }

}

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  int result;
  {
  	PetscInitialize(&argc, &argv, nullptr, nullptr);

  	result = Catch::Session().run(argc, argv);

  	PetscFinalize();
  }
  	Kokkos::finalize();
  	return result;
}	
