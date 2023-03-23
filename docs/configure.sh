#!/bin/zsh
export FFLAGS=" -arch arm64"
export CFLAGS=" -arch arm64"
export LDFLAGS=" -arch arm64"
export CXXFLAGS=" -arch arm64"

./configure \
        --with-fortran-lib="-L/usr/local/gfortran/lib -lgfortran" \
	--without-Love --without-kml --without-Sealevelchange \
	--prefix=$ISSM_DIR \
	--without-wrappers \
	--enable-debugging \
	--enable-development \
	--with-mpi-include="$ISSM_DIR/externalpackages/petsc/install/" \
	--with-mpi-libflags="-L$ISSM_DIR/externalpackages/petsc/install/ -lmpich" \
	--with-petsc-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-metis-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-scalapack-dir="$ISSM_DIR/externalpackages/petsc/install/" \
	--with-mumps-dir="$ISSM_DIR/externalpackages/petsc/install/" \
	--with-ocean="yes" \
	--with-numthreads=4

	exit 1
