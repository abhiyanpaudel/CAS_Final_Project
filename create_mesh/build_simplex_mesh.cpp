#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_mpi.hpp>

int main(int argc, char** argv) {
  Omega_h::Library lib(&argc, &argv);
  Omega_h::Mesh mesh = Omega_h::build_box(
    lib.world(),          // MPI communicator
    OMEGA_H_SIMPLEX,      // Element type: simplex (triangles in 2D)
    1.0, 1.0, 0.0,        // Box dimensions: x=1.0, y=1.0, z=0.0 (2D)
    10, 10, 0,            // Number of elements: nx=10, ny=10, nz=0 (2D)
    false                 // Periodicity: false (non-periodic)
  );

  Omega_h::write_mesh("simplex_mesh", &mesh); // Save mesh to files
  return 0;
}
