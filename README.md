# 🔄 Particle to Mesh Mapping (particle2mesh_map)

![Badge Status: Alpha](https://img.shields.io/badge/Status-Alpha-orange.svg)
![Badge License: BSD 3-Clause](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)

A high-performance parallel utility for mapping particle-based data to finite element meshes using PETSc and Kokkos acceleration.

> [!IMPORTANT]  
> This software is currently in **alpha status**. APIs may change in future releases as functionality is expanded.

## 🔍 Overview

This software implements efficient algorithms for projecting data from distributed particle sets onto finite element meshes using mass-conservative interpolation. It uses GPU acceleration through Kokkos and advanced sparse linear algebra operations via PETSc to enable high-performance mapping between arbitrary particle/mesh distributions.

Key features:
- ⚡ Fast particle-to-mesh (P2M) data interpolation
- 🖥️ GPU acceleration using the Kokkos performance portability library
- 📊 Sparse matrix operations using PETSc's AIJKokkos format
- 🔺 Support for triangular finite element meshes (with plans for tetrahedral and other element types)
- 📐 Barycentric coordinate-based interpolation for accurate field representation

## 📚 Dependencies

Required libraries:
- [Omega_h](https://github.com/sandialabs/omega_h) - Unstructured mesh adaptation library
- [PETSc](https://petsc.org/) - Portable, Extensible Toolkit for Scientific Computation
- [Kokkos](https://github.com/kokkos/kokkos) - Performance portability programming model
- [MeshField](https://github.com/SCOREC/meshFields) - Mesh-based field discretization
- [PCMS](https://github.com/SCOREC/pcms) - Particle to continuum mapping library
- [Gmsh](https://gmsh.info/) - Mesh generation (for creating test meshes)
- MPI - Message Passing Interface

## 🛠️ Installation

### Prerequisites

> [!TIP]
> Make sure to load appropriate modules for your system before installing dependencies. On SCOREC machines, you may need to load gcc, cuda, and mpi modules.

Ensure you have the required dependencies installed. The following paths need to be set:

```bash
# PETSc
export PETSC_DIR=/path/to/petsc/install
export PETSC_ARCH=arch-linux-c-debug

# Omega_h
export OMEGA_H_DIR=/path/to/omega_h/install

# Kokkos
export KOKKOS_DIR=/path/to/kokkos/install

# MeshField
export MESHFIELD_DIR=/path/to/meshfield/install

# PCMS
export PCMS_DIR=/path/to/pcms/install
```

> [!INFO]
> If you're working on a system with pre-installed libraries, check with your system administrator for the correct paths to use.

### Building from Source

```bash
# Clone the repository
git clone https://github.com/username/particle2mesh_map.git
cd particle2mesh_map

# Create build directory
mkdir -p build
cd build

# Configure with CMake
cmake ..

# Build
make
```

> [!TIP]
> Use `make -j N` where N is the number of cores to speed up compilation.

## 🚀 Usage

### Running the Mapping Application

```bash
cd build
./point2MeshMapping <source_mesh_path> <target_mesh_path> [options]
```

#### Required Arguments

- `source_mesh_path`: Path to the source mesh (.osh format)
- `target_mesh_path`: Path to the target mesh (.osh format)

#### Optional Arguments

- `-mat_type aijkokkos`: Use AIJKokkos matrix format (recommended for GPU acceleration)
- `-use_gpu_aware_mpi 0/1`: Enable/disable GPU-aware MPI
- `-ksp_type [solver]`: Specify KSP solver (e.g., cg, gmres)
- `-pc_type [preconditioner]`: Specify preconditioner (e.g., jacobi, bjacobi)

> [!TIP]
> For optimal performance on GPUs, use `-mat_type aijkokkos` and enable GPU-aware MPI if supported by your system.

### Example Workflow

1. Generate meshes using the provided scripts:
```bash
cd create_mesh
./process_geo_files.sh
```

2. Run the mapping operation:
```bash
cd build
./point2MeshMapping ../create_mesh/source_mesh/source_mesh.osh ../create_mesh/target_mesh/target_mesh.osh -mat_type aijkokkos -use_gpu_aware_mpi 0
```

3. Visualization:
The application outputs solution fields in VTK format that can be viewed with visualization tools like ParaView.

> [!CAUTION]
> Ensure your mesh files are properly formatted. Invalid mesh files can lead to runtime errors or incorrect results.

## 🧠 Implementation Details

### Core Components

- `main.cpp`: Entry point with mesh loading and function evaluation
- `massPhiMatrixSolver.hpp`: High-level driver for solving the mass-phi system
- `calculateMassMatrix.hpp`: Constructs the mass matrix for the mesh
- `calculatePhiMatrix.hpp`: Builds the phi (interpolation) matrix
- `massMatrixIntegrator.hpp`: Integrates shape functions for mass matrix computation

> [!INFO]
> The implementation currently focuses on 2D triangular meshes. Future versions will extend support to 3D elements.

### Mathematical Framework

The mapping operation solves the system:
```
M x = Φ f
```
where:
- M is the mass matrix on the target mesh
- Φ is the interpolation matrix from source to target
- f is the value at source points
- x is the resulting field on the target mesh

> [!TIP]
> For large problems, consider using iterative solvers like conjugate gradient (`-ksp_type cg`) with appropriate preconditioners.

## 🧪 Testing and Validation

The provided test meshes in `create_mesh/` directory can be used for verifying the implementation. Running the mapper with these meshes serves as a basic functional test.

> [!INFO]
> More comprehensive test cases will be added in future updates.

## 👥 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

> [!TIP]
> When developing, use the debug print statements in the code for troubleshooting. Add `#define DEBUG` to enable additional output.

## 📝 License

This project is licensed under the BSD 3-Clause License:

```
Copyright (c) 2023-2025, The particle2mesh_map Authors.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## Acknowledgments

- The Omega_h and MeshField teams for their robust mesh handling libraries
- PETSc development team for their high-performance linear algebra framework
- Kokkos team for their performance portability layer
- The developers of the mesh generating scripts in the create_mesh directory that facilitate testing and validation

> [!TIP]
> Check the documentation of each dependency for more detailed information on their usage and capabilities.
