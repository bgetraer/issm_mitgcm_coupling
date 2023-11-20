#!/bin/bash
#INPUT: Models/org
#OUTPUT: Solves and saves to Models/Steadystate
#USEAGE: queue_runsteadystate.sh Models/org

#Set model directory
RUN_DIR=$1
#Set directory path for runcouple script
RUNSTEADYSTATE_DIR="$JPL_DIR/proj-getraer/issm_mitgcm_coupling/bgetraer_pfe/runsteadystate"
#set .queue filename
prefix=$(echo $RUN_DIR | sed -e "s/\// /g" | awk {'print $(NF-1)'})
queue_filename="$RUN_DIR/${prefix}_runsteadystate.queue"

echo "*   - writing .queue file   $queue_filename"

#cd $JPL_DIR/proj-getraer/issm_mitgcm_coupling/bgetraer_pfe/runsteadystate
#cat <<EOF > runsteadystate.queue
cat <<EOF > $queue_filename
#PBS -S /bin/bash
#PBS -l select=1:ncpus=28:model=bro
#PBS -q debug
#PBS -l walltime=0:20:00
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

#Export some variables
export PATH='$PATH:.'
export MPI_LAUNCH_TIMEOUT=800
export MPI_GROUP_MAX=800

#ISSM stuff
export ISSM_DIR='/nobackup/bgetraer/trunk-jpl'
source $ISSM_DIR/etc/environment.sh

#move to the run directory, link the MCC files
cd $RUN_DIR
ln -s $RUNSTEADYSTATE_DIR/mccfiles/run_MCCexecutable.sh ./
ln -s $RUNSTEADYSTATE_DIR/mccfiles/MCCexecutable ./

./run_MCCexecutable.sh $ISSM_DIR/lib:$PETSC_DIR/lib:$MPI_ROOT/lib:${MKLROOT}/lib/intel64_lin:${MKLROOT}/../compiler/lib/intel64_lin:${ISSM_DIR}/externalpackages/triangle/install/lib:/nasa/matlab/2022b $2
EOF


echo 'Sending job to debug queue'
qsub $queue_filename
