#!/bin/bash
#USEAGE: longqueue_runcouple.sh rundir envfile

#Set directory paths for inputs
RUN_DIR=$1
ENVFILE_PATH=$2
#Set directory path for runcouple script and start there
RUNCOUPLE_DIR="$JPL_DIR/proj-getraer/issm_mitgcm_coupling/bgetraer_pfe/runcouple"
#set .queue filename
prefix=$(echo $RUN_DIR | sed -e "s/\// /g" | awk {'print $(NF-1)'})
queue_filename="$RUN_DIR/${prefix}_pickupcouple.queue"

echo "*   - writing .queue file   $queue_filename"
#write the .queue file
cat <<EOF > $queue_filename
#PBS -S /bin/bash
#PBS -l select=1:ncpus=28:mpiprocs=28:model=bro
#PBS -q long
#PBS -l walltime=120:00:00
#PBS -m e
#PBS -W group_list=s2541
#PBS -o $RUN_DIR/run.outlog
#PBS -e $RUN_DIR/run.errlog

. /usr/share/modules/init/bash

#load modules
module load matlab/2022b
module load mpi-hpe/mpt
module load comp-intel
module load petsc/3.17.3_intel_mpt_py
module load hdf5/1.8.18_mpt hdf4/4.2.12 netcdf/4.4.1.1_mpt

#Export some variables
export PATH='$PATH:.'
export MPI_LAUNCH_TIMEOUT=800
export MPI_GROUP_MAX=800

#ISSM stuff
export ISSM_DIR="/nobackup/bgetraer/trunk-jpl"
source $ISSM_DIR/etc/environment.sh

#move to the run directory, link the MCC files
cd $RUN_DIR
ln -s $RUNCOUPLE_DIR/mccfiles/run_MCCexecutable.sh ./
ln -s $RUNCOUPLE_DIR/mccfiles/MCCexecutable ./

#run the runcouple executable with the envfile input
./run_MCCexecutable.sh /nasa/netcdf/4.4.1.1_mpt/lib:$ISSM_DIR/lib:$PETSC_DIR/lib:$MPI_ROOT/lib:${MKLROOT}/lib/intel64_lin:${MKLROOT}/../compiler/lib/intel64_lin:${ISSM_DIR}/externalpackages/triangle/install/lib:/nasa/matlab/2022b $ENVFILE_PATH
EOF

echo '*   - sending job to long queue'
qsub $queue_filename
echo '************************************************************************************'
