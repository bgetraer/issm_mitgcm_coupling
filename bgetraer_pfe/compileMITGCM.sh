#!/bin/bash
#compile.sh is a script that compiles and links MITgcm
#   Compiles a new instance of MITgcm into build_dir using the configuration files in code_dir
#
#   Useage:
#      bash compile.sh model_dir    #compiles from the sub-directories in model_dir
#      bash compile.sh              #compiles from the current sub-directories

#define the model directory
if [ -z "$1" ]   #no input
then
	model_dir=.
else             #read input
	model_dir=$1
fi

#check that sub-directories exist
dirs=(build code)
for dir in "${dirs[@]}"
do
	if  [ ! -d "$model_dir/$dir" ]
	then
		echo "$BASH_SOURCE: line $LINENO: '$model_dir/$dir' is not a directory";
		exit 1
	fi
done

#system directories
genmake2='$ROOTDIR/tools/genmake2'
optfile_dir='$ROOTDIR/tools/build_options/linux_amd64_ifort+mpi_ice_nas'

#compile MITgcm
cd $model_dir/build
$genmake2 -mpi -mo ../code -optfile $optfile_dir -rd $ROOTDIR #generate Makefile
make CLEAN    #remove previous builds in build_dir
make depend   #create symbolic links from the local directory to the source file locations
make          #compile code and create executable file mitgcmuv
