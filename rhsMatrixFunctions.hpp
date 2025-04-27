#ifndef COMPUTING_AT_SCALE_RHS_MATRIX_FUNCTIONS_HPP
#define COMPUTING_AT_SCALE_RHS_MATRIX_FUNCTIONS_HPP

#include <pcms/point_search.h>
#include <Omega_h_shape.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_for.hpp>
#include <MeshField_Shape.hpp>
#include <petscmat.h>
#include <petscvec_kokkos.hpp>
#include <Omega_h_array.hpp>
#include <Kokkos_Core.hpp>


Omega_h::LOs find_element_from_source_points(const Omega_h::Mesh& mesh, const Omega_h::Reals& source_coordinates){
	
	const auto dim = mesh.dim();
	const auto& nfaces = mesh.nfaces();
	const auto& nsources = source_coordinates.size()/dim; 
	
	Kokkos::View<pcms::Real*[2]> source_points("points to find the elements", nsources);
	
	Omega_h::parallel_for(nsources, 
			OMEGA_H_LAMBDA(const int i) {
		
		source_points(i, 0) = source_coordinates[i * dim];
		source_points(i, 1) = source_coordinates[i * dim + 1];

	});

	Kokkos::fence();

	Omega_h::Write<Omega_h::LO> element_ids(nsources, 0, "stores the element ids for corresponding source points"); 
	pcms::GridPointSearch search_element(mesh, 10, 10);
	auto results = search_element(source_points);

	Omega_h::parallel_for(nsources, OMEGA_H_LAMBDA(const int id){
		element_ids[id] = results(id).tri_id;
	});

	return Omega_h::read(element_ids);
}

struct Results{
	Kokkos::View<PetscInt*[3]> row_ids;
	Kokkos::View<PetscInt*[3]> vals;
};

Results get_row_ids_and_values(const Omega_h::Mesh& mesh, const Omega_h::Reals& source_coordinates){
	
	const auto dim = mesh.dim(); // 2
	const auto& nfaces = mesh.nfaces();
	const auto& nsources = source_coordinates.size()/dim;
	const auto& faces2nodes = mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;
	const auto& coordinates = mesh.coords();

	// get element ids for corrponding source points 
	auto element_ids = find_element_from_source_points(mesh, source_coordinates);

	// allocate per-column metadata 
	Kokkos::View<int*[3]> row_ids("row indices", nsources, 3); 
	Kokkos::View<int*[3]> vals("nonzero values", nsources, 3);
	
	// fill row_ids and vals
	Omega_h::parallel_for(nsources, OMEGA_H_LAMBDA(const int id) {	
		
		auto current_src_el_id = element_ids[id];
		const auto current_el_verts = Omega_h::gather_verts<3>(faces2nodes, current_src_el_id);
	
		// row indices = 3 vertex ids
		for (int i = 0; i < 3; ++i){
			row_ids(id, i) = current_el_verts[i];			
		}

		const Omega_h::Few<Omega_h::Vector<2>,3> current_el_vert_coords = 
				Omega_h::gather_vectors<3,2>(coordinates, current_el_verts);
		Omega_h::Vector<2> global_coordinate{source_coordinates[id * dim + 0], source_coordinates[id * dim + 1]};
		auto barycentric_coordinate = barycentric_from_global<2,3>(global_coordinate, current_el_vert_coords);		
		
		Kokkos::Array<Omega_h::Real, 3> barycentric_coords_array;
		
		for (int j = 0; j < 3; ++j) {
			barycentric_coords_array[j] = barycentric_coordinate[j];
		}

		MeshField::LinearTriangleShape lts;
		auto phi = lts.getValues(barycentric_coords_array);

		for (int k = 0; k < 3; ++k) {
			vals(id, k) = phi[k];	
		}

	});

	Kokkos::fence();
	return Results{row_ids, vals); 

}

Mat createCSRforRHSMatrix(const Omega_h::Mesh& mesh, const Omega_h::Reals& source_coordinates){
	// get per-column data
	auto res = get_row_ids_and_values(mesh, source_coordinates);
	PetscInt nrows = mesh.nverts();
	PetscInt ncols = source_coordinates.size()/mesh.dim();
	PetscInt nnz = 3 * ncols;

	// build ia
	Kokkos::View<PetscInt*> ia("ia", nrows+1);
	Kokkos::deep_copy(ia_d, 0);
	Omega_h::parallel_for(ncols, OMEGA_H_LAMBDA(const int i) {
		for (int j = 0; j < 3; ++j) {
			Kokkos::atomic_increment(&ia(res.row_ids(i,j) + 1));	
		}
	});
	Kokkos::fence();

	int result = 0;
	Kokkos::parallel_scan("inclusive sum", nrows+1, KOKKOS_LAMBDA(const int i, PetscInt& inner_sum, bool final){
		inner_sum += ia(i);
		if (final){
			ia(i) = inner_sum;
		}
	},result);

	Kokkos::fence();

	// allocate ja, a 
	Kokkos::View<PetscInt*> ja("ja", nnz);
	Kokkos::View<PetscScalar*> a("ja", nnz);
	Kokkos::View<PetscInt*> next("next", nrows);
	
	// init next
	Omega_h::parallel_for(nrows, OMEGA_H_LAMBDA(int i) {
		next(i) = ia(i);
	});
	Kokkos::fence();

	// scatter 
	Omega_h::parallel_for(ncols, OMEGA_H_LAMBDA(const int j) {
		for (int k = 0; k < 3; ++k){
			auto r = res.row_ids(j,k);
			auto slot = Kokkos::atomic_fetch_add(&next(r), 1);
			ja(slot) = j;
			a(slot) = res.vals(j,k);
		}
	});
	Kokkos::fence();

	// PETSc CSR matrix
	Mat A;
	PetscCall(MatCreateSeqAIJKokkosWithKokkosViews(
				PETSC_COMM_SELF,
				nrows,
				ncols,
				ia_d, ja_d, a_d,
				&A));
	return A;
}











#endif 
