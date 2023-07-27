#!/bin/bash

#system directories
genmake2=$ROOTDIR/tools/genmake2
optfile_dir=$ROOTDIR/tools/build_options/linux_amd64_gfortran

#compile MITgcm
cd ./build
rm *
$genmake2 -mpi -mo ../code -optfile $optfile_dir  -rd $ROOTDIR #generate Makefile
make CLEAN    #remove previous builds in build_dir
make depend   #create symbolic links from the local directory to the source file locations
make          #compile code and create executable file mitgcmuv
cd ..
