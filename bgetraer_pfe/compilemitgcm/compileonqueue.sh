#!/bin/bash
cat <<EOF > compile.queue
#PBS -S /bin/bash
#PBS -l select=1:ncpus=20:model=bro
#PBS -q devel
#PBS -l walltime=0:20:00
#PBS -m e
#PBS -W group_list=s2541
#PBS -o /nobackup/bgetraer/issmjpl/proj-getraer/run.outlog
#PBS -e /nobackup/bgetraer/issmjpl/proj-getraer/run.errlog

#. /usr/share/modules/init/bash
#
##load modules
module load comp-intel/2020.4.304
module load /nasa/intel/impi/2021.3/modulefiles/mpi/2021.3.0
#module load matlab/2022b

#MITgcm stuff
#umask 027
#limit stacksize unlimited

#Export some variables
#export PATH='$PATH:.'
#export MPI_LAUNCH_TIMEOUT=800
#export MPI_GROUP_MAX=800

#ISSM stuff
#export ISSM_DIR='/nobackup/bgetraer/trunk-jpl'
#source $ISSM_DIR/etc/environment.sh

#move and start simulation
source ~/.bashrc
./compile.sh ../
EOF

echo 'Sending job to the queue'
qsub compile.queue
