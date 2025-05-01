source /lore/paudea/scripts/loads-rhel9.sh


export DEVICE_ARCH=ADA89-meshField


SOURCE_DIR="${SOURCE_DIR:-/lore/paudea/sources}"

BUILD_DIR="${BUILD_DIR:-/lore/paudea/build}"


DEVICE_ARCH="${DEVICE_ARCH:-ADA89-meshField}"

DEPENDENCY_DIR="${DEPENDENCY_DIR:-/lore/mersoj2/laces-software/build}"

CURDIR=$PWD

export NVCC_WRAPPER_DEFAULT_COMPILER=`which mpicxx`

#-DCMAKE_CXX_COMPILER=$BUILD_DIR/${DEVICE_ARCH}/kokkos/install/bin/nvcc_wrapper \

#	-DCMAKE_CXX_COMPILER=`which mpicxx` \
cmake -S . -B build \
	-DCMAKE_VERBOSE_MAKEFILE=ON \
	-DCMAKE_CXX_COMPILER=$BUILD_DIR/${DEVICE_ARCH}/kokkos-meshField/install/bin/nvcc_wrapper \
    -DCMAKE_C_COMPILER=`which mpicc` \
    -DOmega_h_USE_Kokkos=ON \
    -DOmega_h_USE_CUDA=ON \
    -DOmega_h_DIR=$BUILD_DIR/${DEVICE_ARCH}/omega_h-meshField/install/lib64/cmake/Omega_h/ \
    -DKokkos_ROOT=$BUILD_DIR/${DEVICE_ARCH}/kokkos-meshField/install/ \
	-Dpcms_ROOT=$BUILD_DIR/${DEVICE_ARCH}/pcms-meshField/install/ \
    -Dperfstubs_DIR=$DEPENDENCY_DIR/perfstubs/install/lib/cmake/ \
    -DADIOS2_ROOT=$DEPENDENCY_DIR/adios2/install/ \
    -DCMAKE_BUILD_TYPE=Debug \
	-Dmeshfields_DIR=$BUILD_DIR/${DEVICE_ARCH}/meshField/install/lib64/cmake/meshfields/ \
	-DPETSC_ARCH=$SOURCE_DIR/petsc/linux-gnu-gpu-kokkos/ \
	-DPETSC_DIR=$SOURCE_DIR/petsc/ \
    -DCatch2_ROOT=$DEPENDENCY_DIR/Catch2/install/ 

cmake --build build -j8
cd $CURDIR
