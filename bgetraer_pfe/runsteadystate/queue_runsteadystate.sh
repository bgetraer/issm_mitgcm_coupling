#!/bin/bash
#INPUT: Models/org
#OUTPUT: Solves and saves to Models/Steadystate
#USEAGE: queue_runsteadystate.sh Models/org

cd $JPL_DIR/proj-getraer/issm_mitgcm_coupling/bgetraer_pfe/runsteadystate

cat <<EOF > runsteadystate.queue
#PBS -S /bin/bash
#PBS -l select=1:ncpus=28:model=bro
#PBS -q devel
#PBS -l walltime=0:20:00
#PBS -m e
#PBS -W group_list=s2541
#PBS -o ./run.outlog
#PBS -e ./run.errlog

. /usr/share/modules/init/bash

#load modules
module load matlab/2022b
module load mpi-hpe/mpt
module load comp-intel
module load petsc/3.17.3_intel_mpt_py

#Export some variables
export PATH='$PATH:.'
export MPI_LAUNCH_TIMEOUT=800
export MPI_GROUP_MAX=800

#ISSM stuff
export ISSM_DIR='/nobackup/bgetraer/trunk-jpl'
source $ISSM_DIR/etc/environment.sh

./mccfiles/run_MCCexecutable.sh $ISSM_DIR/lib:$PETSC_DIR/lib:$MPI_ROOT/lib:${MKLROOT}/lib/intel64_lin:${MKLROOT}/../compiler/lib/intel64_lin:${ISSM_DIR}/externalpackages/triangle/install/lib:/nasa/matlab/2022b $1
EOF


echo 'Sending job to the queue'
qsub runsteadystate.queue
