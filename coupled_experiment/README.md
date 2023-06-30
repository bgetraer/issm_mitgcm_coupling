### To get a copy of this repo
  git clone git@github.com:MITgcm-contrib/issm_mitgcm_coupling
  cd issm_mitgcm_coupling
  export BASEDIR=`pwd`
  export EXPDIR=$BASEDIR/coupled_experiment

### To get and/or build ISSM
# Follow instructions in https://issm.jpl.nasa.gov/download
# Instructions for installing on mac silicon are in
# $BASEDIR/issm_mitgcm_coupling/doc/issm_on_mac_silicon.txt

### To build MITgcm
# To make the experiment interesting, we clone Dan's github repo,
# checkout `branch_runoff_2`, which is based on checkpoint 67v and
# includes parameters in shelfice relating to subglacial runoff,
# set the ROOTDIR environment variable, and build MITgcm.
  git clone git@github.com:dngoldberg/MITgcm
  cd MITgcm
  git checkout branch_runoff_2
  export ROOTDIR=`pwd`

# To compile on linux cluster
  cd $EXPDIR/build
  $ROOTDIR/tools/genmake2 -mods=../code -mpi -of=$ROOTDIR/tools/build_options/linux_amd64_gfortran --rootdir=$ROOTDIR

# To compile on mac
  cd $EXPDIR/build
  $ROOTDIR/tools/genmake2 -mods=../code -mpi -of=$ROOTDIR/tools/build_options/darwin_amd64_gfortran --rootdir=$ROOTDIR

# Compile and link code
  make depend
  make -j

### To run ISSM/MITgcm
  cd $EXPDIR
  matlab
  addpath input

# To run on mac with np=8
  nPx=2; nPy=4;

# To prepare the initial ice-shelf state
  steps=1:3; runme

# To list available experiments
# step # 1: 'GetMITgcm'
# step # 2: 'Mesh_generation'
# step # 3: 'SteadystateNoSlip'
# step # 4: 'RunCouple'
# step # 5: 'RunCouple2'
# step # 6: 'RunCouple3'
# step # 7: 'TestDrift'
# step # 8: 'TestDrift2'
  org

# To run Dan's matlab-based experiment
  steps=4; runme
 
# To run Helene's matlab-based experiment
  steps=5; runme
 
# To run MPI-based experiment
  steps=6; runme

# To look at output
  steps=8; runme

# Note that there are self-explanatory run time parameters at the
# top of the `'RunSingleCoupleStep'` section.

### To test the correction mechanism
# `test.m` will run the coupled step (step 4) twice, one with `alpha_correction`
#   set to 0 and then with it set to 1. (It will save the run directories and
#   model object from each run.) *need to do steps 1 thru 3 first*.
# `test2.m` will then do the testing step (step 5) for each test. Currently,
#   only the last coupled time step is analysed in `runme.m` but it could be
#   adapted to show a time series of issm-shelfice mass misfit. Note the `nsteps`
#   variable in `test2.m` must be the same as `test.m`.
