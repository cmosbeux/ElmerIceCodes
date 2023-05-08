#!/bin/bash

nb_params=$#
if [ $nb_params -lt 4 ];then
        echo "USAGE sh $0 MPILIB(openmpi/intelmpi)   COMPILER (gnu/intel) COMPILER_VERSION  CLUSTER (fc18/nautilus/centos7/lachouf/debian/froggy) "
        echo "Example: bash  ige_install_elmerfem_new.sh  intelmpi intel 18 dahu 71c874d"
        exit
fi

set -x 

MPILIB=$1
COMPILER=$2
COMPILER_VERSION=$3
CLUSTER=$4
ELMERVERSION=$5

INSTALL_DIR=/home/mchekki/LibGlace/versions

INSTALL_DIR_BETTIK=/bettik/mosbeuxc/INSTALL/ELMER/

ELMERFEM_VERSION=elmerfemdev-${ELMERVERSION}
MUMPS_VERSION=mumps-5.2.0
MUMPS_VERSION=mumps-5.2.0-mkl
SCALAPACK_VERSION=scalapack-2.0.2
LAPACK_VERSION=lapack-3.8.0
SCOTCH_VERSION=scotch-6.0.7
METIS_VERSION=metis-5.1.0
NETCDF_VERSION=netcdf-4.7.2
HDF5_VERSION=hdf5-1.10.4
XIOS_VERSION=xios-2.0_rev2320  #xios-2.0_rev1689
ZLIB_VERSION=zlib-1.2.11
SZIP_VERSION=szip-2.1.1
MMG_VERSION=mmg-5.6.0
NNC_VERSION=nnc-1.85
CSAC_VERSION=csac-1.22
FORTRANGIS_VERSION=fortrangis-2.5
PROJ_VERSION=proj-4.9.3
HYPRE_VERSION=hypre-2.1


#
if [ ! -d ${INSTALL_DIR} ]; then

echo "Folder ${INSTALL_DIR} does not exit"
echo "Make sure you have chosen an existing cluster"

exit
fi

# Global Environment variables

export FC=mpiifort
export CC=mpiicc
export CXX=mpiicpc


#Environment variable


export ELMERFEM_DIR=${INSTALL_DIR_BETTIK}/elmerfem/${ELMERFEM_VERSION}"_"${COMPILER}${COMPILER_VERSION}
export ELMERFEM_LIB=${ELMERFEM_DIR}/lib

export MUMPS_DIR=${INSTALL_DIR}/mumps/${MUMPS_VERSION}"_"${COMPILER}${COMPILER_VERSION}
export MUMPS_INC=${MUMPS_DIR}/include
export MUMPS_LIB=${MUMPS_DIR}/lib

export SCALAPACK_DIR=${INSTALL_DIR}/scalapack/${SCALAPACK_VERSION}"_"${COMPILER}${COMPILER_VERSION}
export SCALAPACK_DIR=/home/mchekki/.nix-profile/mkl
export SCALAPACK_INC=${SCALAPACK_DIR}/include
export SCALAPACK_LIB=${SCALAPACK_DIR}/lib/intel64


export METIS_DIR=${INSTALL_DIR}/metis/${METIS_VERSION}"_"${COMPILER}${COMPILER_VERSION}
export METIS_INC=${METIS_DIR}/include
export METIS_LIB=${METIS_DIR}/lib
export metis_libs=""
export metis_libs="-L${METIS_LIB} -lparmetis -lmetis"


export LAPACK_DIR=${INSTALL_DIR}/lapack/${LAPACK_VERSION}"_"${COMPILER}${COMPILER_VERSION}
export LAPACK_DIR=/home/mchekki/.nix-profile/mkl
export LAPACK_INC=${LAPACK_DIR}/include
export LAPACK_LIB=${LAPACK_DIR}/lib

export SCOTCH_DIR=${INSTALL_DIR}/scotch/${SCOTCH_VERSION}"_"${COMPILER}${COMPILER_VERSION}
export SCOTCH_INC=${SCOTCH_DIR}/include
export SCOTCH_LIB=${SCOTCH_DIR}/lib
export scotch_libs=""
export scotch_libs="-L${SCOTCH_LIB} -lptscotcherr -lptscotcherrexit  -lptesmumps -lptscotch -lscotcherr  -lscotcherrexit  -lscotch"

export NETCDF_DIR=${INSTALL_DIR}/netcdf/${MPILIB}/${NETCDF_VERSION}"_"${COMPILER}${COMPILER_VERSION}
export NETCDFROOT=${NETCDF_DIR}
export NETCDF_INC=${NETCDF_DIR}/include
export NETCDF_LIB=${NETCDF_DIR}/lib

export HDF5_DIR=${INSTALL_DIR}/hdf5/${MPILIB}/${HDF5_VERSION}"_"${COMPILER}${COMPILER_VERSION}
export HDF5ROOT=${HDF5_DIR}
export HDF5_INC=${HDF5_DIR}/include
export HDF5_LIB=${HDF5_DIR}/lib

export XIOS_DIR=${INSTALL_DIR}/xios/${MPILIB}/${XIOS_VERSION}"_"${COMPILER}${COMPILER_VERSION}
export XIOSROOT=${XIOS_DIR}
export XIOS_INC=${XIOS_DIR}/include
export XIOS_LIB=${XIOS_DIR}/lib

export ZLIB_DIR=${INSTALL_DIR}/zlib/${ZLIB_VERSION}"_"${COMPILER}${COMPILER_VERSION}
export ZLIB_LIB=${ZLIB_DIR}/lib

export SZIP_DIR=${INSTALL_DIR}/szip/${SZIP_VERSION}"_"${COMPILER}${COMPILER_VERSION}
export SZIP_LIB=${SZIP_DIR}/lib

export MMG_DIR=${INSTALL_DIR}/mmg/${MMG_VERSION}"_"${COMPILER}${COMPILER_VERSION}
export MMGROOT=${MMG_DIR}
export MMG_INC=${MMG_DIR}/include
export MMG_LIB=${MMG_DIR}/lib

export NNC_DIR=${INSTALL_DIR}/nnc/${NNC_VERSION}"_"${COMPILER}${COMPILER_VERSION}
export NNC_INC=${NNC_DIR}/include
export NNC_LIB=${NNC_DIR}/lib

export CSAC_DIR=${INSTALL_DIR}/csac/${CSAC_VERSION}"_"${COMPILER}${COMPILER_VERSION}
export CSAC_INC=${CSAC_DIR}/include
export CSAC_LIB=${CSAC_DIR}/lib

export FORTRANGIS_DIR=${INSTALL_DIR}/fortrangis/${FORTRANGIS_VERSION}"_"${COMPILER}${COMPILER_VERSION}
export FORTRANGISROOT=${FORTRANGIS_DIR}
export FORTRANGIS_INC=${FORTRANGIS_DIR}/include
export FORTRANGIS_LIB=${FORTRANGIS_DIR}/lib

export PROJ_DIR=${INSTALL_DIR}/proj/${PROJ_VERSION}"_"${COMPILER}${COMPILER_VERSION}
export PROJROOT=${PROJ_DIR}
export PROJ_INC=${PROJ_DIR}/include
export PROJ_LIB=${PROJ_DIR}/lib

export HYPRE_DIR=${INSTALL_DIR}/hypre/${HYPRE_VERSION}"_"${COMPILER}${COMPILER_VERSION}
export HYPREROOT=${HYPRE_DIR}
export HYPRE_INC=${HYPRE_DIR}/include
export HYPRE_LIB=${HYPRE_DIR}/lib


export LD_LIBRARY_PATH=\
${MUMPS_LIB}:\
${SCALAPACK_LIB}:\
${LAPACK_LIB}:\
${NETCDF_LIB}:\
${HDF5_LIB}:\
${ZLIB_LIB}:\
${SZIP_LIB}:\
${MMG_LIB}:\
${NNC_LIB}:\
${CSAC_LIB}:\
${HYPRE_LIB}:\
${ELMERFEM_LIB}/elmersolver:\
$LD_LIBRARY_PATH


# Example build-script for Elmer on Linux

export LANG=C
ulimit -s unlimited
ulimit -v unlimited

CMAKE=cmake


# BUILD ELMER 
TIMESTAMP=$(date +"%m-%d-%y")
ELMERFEM_SRC=/bettik/mosbeuxc/INSTALL/ELMER/Elmer_${ELMERVERSION}_SRC
BUILD_DIR="$ELMERFEM_SRC/builddir"

rm -rf $BUILD_DIR ; mkdir -p $BUILD_DIR

echo "Building Elmer from within " ${BUILD_DIR}
echo "installation into " ${ELMERFEM_DIR}
cd ${BUILD_DIR}
pwd
ls -ltr

 #-DCMAKE_C_FLAGS="-O3 -xHost -ip -xAVX2 -qopenmp  -traceback" \
 #-DCMAKE_Fortran_FLAGS="-O3 -qopenmp  -traceback" \
$CMAKE  -Wno-dev  $ELMERFEM_SRC -DCMAKE_INSTALL_PREFIX=$ELMERFEM_DIR \
 -DCMAKE_Fortran_COMPILER=$FC \
 -DCMAKE_C_COMPILER=$CC \
 -DCMAKE_BUILD_TYPE="Release" \
 -DCMAKE_C_FLAGS="-O3  -qopenmp  -traceback " \
 -DCMAKE_CXX_COMPILER=$CXX \
 -DCMAKE_Fortran_LINKER=$FC \
 -DCMAKE_Fortran_FLAGS="-cpp -fp-model=consistent -O2 -qopenmp  -traceback " \
 -DWITH_LUA:BOOL=TRUE -DWITH_MKL:BOOL=TRUE \
 -DWITH_OpenMP:BOOL=TRUE \
 -DWITH_ElmerIce:BOOL=TRUE \
 -DWITH_XIOS:BOOL=FALSE \
 -DWITH_MPI:BOOL=TRUE \
 -DMPI_C_COMPILER=$CC -DMPI_CXX_COMPILER=$CXX -DMPI_Fortran_COMPILER=$FC \
 -DWITH_Mumps:BOOL=TRUE \
 -DMumps_LIBRARIES="-L${MUMPS_LIB} -ldmumps -lmumps_common -lpord ${scotch_libs} ${metis_libs} -L${SCALAPACK_LIB} -lmkl_scalapack_lp64   -lmkl_blacs_intelmpi_lp64 -L${ZLIB_LIB} -lz -L${SZIP_LIB} -lsz"\
 -DMumps_INCLUDE_DIR=${MUMPS_INC} \
 -DWITH_ScatteredDataInterpolator:BOOL=TRUE \
 -DWITH_NETCDF:BOOL=TRUE \
 -DNetCDF_INCLUDE_DIR=${NETCDF_INC} -DNetCDF_DIR=${NETCDF_DIR}  \
 -DNetCDF_LIBRARY="${NETCDF_LIB}/libnetcdf.a" \
 -DNetCDFF_LIBRARY="${NETCDF_LIB}/libnetcdff.a"\
 -DPHDF5_INCLUDE_DIR=${HDF5_INC} -DPHDF5_DIR=${HDF5_DIR}  \
 -DPHDF5_LIBRARY="${HDF5_LIB}/libhdf5.a" \
 -DPHDF5HL_LIBRARY="${HDF5_LIB}/libhdf5_hl.a" \
 -DXIOS_INCLUDE_DIR=${XIOS_INC} -DXIOS_DIR=${XIOS_DIR}  \
 -DXIOS_LIBRARY="${XIOS_LIB}/libxios.a" \
 -DWITH_Hypre:BOOL=TRUE \
 -DHYPRE_INCLUDE_DIR=${HYPRE_INC} \
 -DHYPRE_LIBRARY="${HYPRE_LIB}/libHYPRE.a" \
 -DFORTRANGIS_INCLUDE_DIR=${FORTRANGIS_INC} -DFORTRANGIS_DIR=${FORTRANGIS_DIR}  \
 -DFORTRANGIS_LIBRARIES="${FORTRANGIS_LIB}/libfortrangis.a ; ${FORTRANGIS_LIB}/libfortranc.a" \
 -DPROJ_INCLUDE_DIR=${PROJ_INC} -DPROJ_DIR=${PROJ_DIR}  \
 -DPROJ_LIBRARY="${PROJ_LIB}/libproj.a" \
 -DWITH_MMG:BOOL=FALSE \
 -DMMG_INCLUDE_DIR=${MMG_INC} \
 -DMMG_LIBRARY="${MMG_LIB}/libmmg.a" \
 -DCSA_INCLUDE_DIR=${CSAC_INC} \
 -DCSA_LIBRARY="${CSAC_LIB}/libcsa.so" \
 -DNN_INCLUDE_DIR="${NNC_INC}" \
 -DNN_LIBRARY="${NNC_LIB}/libnn.so"  2>&1 | tee cmake.log

make -j 16 && make install  2>&1   |tee make_and_install.log
