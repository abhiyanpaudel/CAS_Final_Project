/**
 * @file calculatePhiMatrix.hpp
 * @brief Functions to calculate Phi matrix 
 * @author [Abhiyan Paaudel]
 * @date May 1, 2025
 *
 * This file contains functions for creating the Phi matrix 
 * The Phi matrix represents the evaluation of basis functions at particle positions.
 */

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
#include <cmath> // For std::abs

/**
 * @brief Finds the mesh elements containing each source point
 * 
 * Uses a grid-based point search algorithm to efficiently locate the 
 * triangle element containing each source point in the mesh.
 * 
 * @param mesh The target Omega_h mesh
 * @param source_coordinates Coordinates of the source points
 * @return Omega_h::LOs Array of element IDs corresponding to each source point
 */
Omega_h::LOs find_element_from_source_points(Omega_h::Mesh& mesh, const Omega_h::Reals& source_coordinates){
	
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

	Omega_h::Write<Omega_h::LO> element_ids(nsources, -1, "stores the element ids for corresponding source points"); 
	pcms::GridPointSearch search_element(mesh, 10, 10);
	auto results = search_element(source_points);

	Omega_h::parallel_for(nsources, OMEGA_H_LAMBDA(const int id){
		element_ids[id] = results(id).tri_id;
	});

	return Omega_h::read(element_ids);
}

/**
 * @brief Structure to hold row indices and values for matrix construction
 * 
 * Contains arrays for the row indices (vertices) and the corresponding 
 * values (shape function evaluations) used in the Phi matrix construction.
 */
struct Results{
	Kokkos::View<PetscInt*[3]> row_ids;    /**< Row indices for each source point and vertex */
	Kokkos::View<PetscScalar*[3]> vals;    /**< Shape function values for each source point and vertex */
};

/**
 * @brief Computes row indices and shape function values for each source point
 * 
 * For each source point, finds the containing element, identifies the element's
 * vertices (row indices), and computes the barycentric coordinates (shape function values)
 * of the source point within that element.
 * 
 * @param mesh The target Omega_h mesh
 * @param source_coordinates Coordinates of the source points
 * @return Results Structure containing row indices and shape function values
 */
Results get_row_ids_and_values(Omega_h::Mesh& mesh, const Omega_h::Reals& source_coordinates){
	
	std::cout << "DEBUG: get_row_ids_and_values - Start" << std::endl;
	const auto dim = mesh.dim(); // 2
	const auto& nfaces = mesh.nfaces();
	const auto& nsources = source_coordinates.size()/dim;
	std::cout << "DEBUG: Number of source points: " << nsources << std::endl;
	
	const auto& faces2nodes = mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;
	const auto& coordinates = mesh.coords();

	// get element ids for corrponding source points 
	std::cout << "DEBUG: Finding elements from source points" << std::endl;
	auto element_ids = find_element_from_source_points(mesh, source_coordinates);
	std::cout << "DEBUG: Element IDs found" << std::endl;
	

	// allocate per-column metadata 
	std::cout << "DEBUG: Allocating row_ids and vals arrays" << std::endl;
	Kokkos::View<PetscInt*[3]> row_ids("row indices", nsources, 3); 
	Kokkos::View<PetscScalar*[3]> vals("nonzero values", nsources, 3);
	
	// fill row_ids and vals 
	std::cout << "DEBUG: Computing shape functions and filling arrays" << std::endl;
	Omega_h::parallel_for(nsources, OMEGA_H_LAMBDA(const int id) {	
		
		auto current_src_el_id = element_ids[id];
		
		const auto current_el_verts = Omega_h::gather_verts<3>(faces2nodes, current_src_el_id);
	
		// row indices = 3 vertex ids
		for (int j = 0; j < 3; ++j){
			row_ids(id, j) = current_el_verts[j];			
		}

		const Omega_h::Few<Omega_h::Vector<2>,3> current_el_vert_coords = 
				Omega_h::gather_vectors<3,2>(coordinates, current_el_verts);
		Omega_h::Vector<2> global_coordinate{source_coordinates[id * dim + 0], source_coordinates[id * dim + 1]};
		
		// Calculate barycentric coordinates
		auto barycentric_coordinate = barycentric_from_global<2,2>(global_coordinate, current_el_vert_coords);
		
		if (id < 10) {  // Only print for the first 10 points to avoid too much output
			printf("Point %d barycentric: [%f, %f, %f]\n", 
				id, 
				barycentric_coordinate[0], 
				barycentric_coordinate[1], 
				barycentric_coordinate[2]);
		}
		
		// Create a Kokkos::Array for the shape function
		//Kokkos::Array<Omega_h::Real, 3> barycentric_coords_array;
		//for (int j = 0; j < 3; ++j) {
	//		barycentric_coords_array[j] = barycentric_coordinate[j];
	//	}

		// Use the shape function
		//MeshField::LinearTriangleShape lts;
		//auto phi = lts.getValues(barycentric_coords_array);
		
		for (int k = 0; k < 3; ++k) {
			vals(id, k) = barycentric_coordinate[k];	
		}
	});

	Kokkos::fence();
	std::cout << "DEBUG: get_row_ids_and_values - Completed" << std::endl;
	return Results{row_ids, vals}; 
}

/**
 * @brief Calculates the Phi matrix for mapping from source points to a target mesh
 *
 * Creates a sparse matrix that is used to map field values from source points to the target mesh.
 * The matrix is constructed in CSR format using row indices and values from shape functions.
 * This implements a finite element projection of point data onto the mesh.
 * 
 * Algorithm:
 * 1. Get row indices and values for each source point
 * 2. Count non-zeros per row and create CSR row pointers (ia array)
 * 3. Fill column indices (ja array) and values (a array)
 * 4. Create the PETSc matrix with these arrays
 * 
 * @param mesh The target Omega_h mesh
 * @param source_coordinates Coordinates of the source points
 * @param[out] Phi_out Pointer to the resulting Phi matrix
 * @return PetscErrorCode PETSc error code (PETSC_SUCCESS if successful)
 */
inline PetscErrorCode calculatePhiMatrix(
    Omega_h::Mesh& mesh,
    const Omega_h::Reals& source_coordinates,
    Mat *Phi_out)
{
  std::cout << "DEBUG: calculatePhiMatrix - Start\n";
  auto res = get_row_ids_and_values(mesh, source_coordinates);
  PetscInt nrows = mesh.nverts();
  PetscInt ncols = source_coordinates.size()/mesh.dim();
  PetscInt nnz   = 3*ncols;
  std::cout << "DEBUG:  nrows="<<nrows<<" ncols="<<ncols<<" nnz="<<nnz<<"\n";

  Kokkos::View<PetscInt*> temp("temp for ia", nrows);
  Kokkos::deep_copy(temp, 0);
  Omega_h::parallel_for(ncols, OMEGA_H_LAMBDA(int i) {
    for (int j = 0; j < 3; ++j) {
      auto r = res.row_ids(i,j);
      if (r < 0 || r >= nrows) {
        printf("ERROR: out-of-range row_id(%d,%d) = %d\n",
               i, j, int(r));
      }
      Kokkos::atomic_increment(&temp(r));
    }
  });
  Kokkos::fence();

  Kokkos::View<PetscInt*> ia("stores scan of temp", nrows+1);
  Kokkos::deep_copy(ia, 0);
  int total_nnz_values = 0;
  Kokkos::parallel_scan("inclusive sum", nrows,
    KOKKOS_LAMBDA(int i, PetscInt &update, bool final) {
      update += temp(i);
      if (final) ia(i + 1) = update;
    }, total_nnz_values);
  Kokkos::fence();
  
  std::cout << "DEBUG: total_nnz_values = " << total_nnz_values
            << "  expected nnz=" << nnz << "\n";
  if (total_nnz_values != nnz) {
    std::cerr << "ERROR: mismatched total nonzeros!\n";
  }
  auto h_ia = Kokkos::create_mirror_view(ia);
  Kokkos::deep_copy(h_ia, ia);
  std::cout << "DEBUG: ia prefix [0..10] = ";
  for (int i = 0; i <= std::min(nrows,10); ++i) std::cout<<h_ia[i]<<" ";
  std::cout<<"\n";

  Kokkos::View<PetscInt*> ja("ja", nnz);
  Kokkos::View<PetscScalar*> a("a", nnz);
  Kokkos::View<PetscInt*>    next("next", nrows);
  Omega_h::parallel_for(nrows, OMEGA_H_LAMBDA(int i) {
    next(i) = ia(i);
  });
  Kokkos::fence();

  Omega_h::parallel_for(ncols, OMEGA_H_LAMBDA(int i) {
    for (int k = 0; k < 3; ++k) {
      auto r = res.row_ids(i,k);
      auto slot = Kokkos::atomic_fetch_add(&next(r), 1);
      ja(slot) = i;
      a(slot)  = res.vals(i,k);
    }
  });
  Kokkos::fence();

  auto h_ja = Kokkos::create_mirror_view(ja);
  auto h_a  = Kokkos::create_mirror_view(a);
  Kokkos::deep_copy(h_ja, ja);
  Kokkos::deep_copy(h_a,  a);
  std::cout << "DEBUG: size of ja :" << h_ja.size() <<"\n";
  std::cout << "DEBUG: size of a :" << h_a.size() <<"\n";
  std::cout << "DEBUG: first 5 (ja, a) = ";
  for (int i = 0; i < std::min(nnz,5); ++i) {
    std::cout << "(" << h_ja[i] << "," << h_a[i] << ") ";
  }
  std::cout << "\n";
  for (int i = 0; i < nnz; ++i) {
    if (h_ja[i] < 0 ) {
       std::cerr << "ERROR: negative entry at i="<<i
                  << " ja="<<h_ja[i]<<" a="<<h_a[i]<<"\n";
        break;
    }
  }

  std::cout<<"DEBUG: creating PETSc matrix\n";
  Mat Phi;
  PetscCall(MatCreateSeqAIJKokkosWithKokkosViews(
    PETSC_COMM_SELF, nrows, ncols, ia, ja, a, &Phi));
  *Phi_out = Phi;
  std::cout<<"DEBUG: calculatePhiMatrix - Completed\n";
  return PETSC_SUCCESS;
}
#endif
