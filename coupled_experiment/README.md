### To get a copy of this repo
  git clone git@github.com:MITgcm-contrib/issm_mitgcm_coupling
  export BASEDIR=`pwd`
  export EXPDIR=$BASEDIR/issm_mitgcm_coupling/coupled_experiment

### To get and/or build ISSM
# Follow instructions here: https://issm.jpl.nasa.gov/download
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
  cd $EXPDIR/build
  $ROOTDIR/tools/genmake2 -mods=../code -mpi -of=$ROOTDIR/tools/build_options/linux_amd64_gfortran --rootdir=$ROOTDIR
# or for compiling on mac:
# $ROOTDIR/tools/genmake2 -mods=../code -mpi -of=$ROOTDIR/tools/build_options/darwin_amd64_gfortran --rootdir=$ROOTDIR
  make depend
  make -j

### To run ISSM/MITgcm open matlab, and run `runme`
  cd $EXPDIR
  matlab
  runme

# To prepare the initial ice-shelf state, 
# change line 5 of `$EXPDIR/runme.m` to `steps=1:3;`

# To run the experiment,
# change line 5 of `$EXPDIR/runme.m` to `steps=4;`

# Note there are self-explanatory run time parameters at the top of the
# `'RunSingleCoupleStep'` section

# On mac or laptop, need to reduce the number of processes (nPx,nPy)
# in SIZE.h before MITgcm compilation and in runme.m
