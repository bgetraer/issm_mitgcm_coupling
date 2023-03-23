### TO BUILD MITgcm

There are parameters in shelfice relating to inputting subglacial runoff. (I suggest we keep these to make the experiment interesting.)

For this reason, please clone my github repo:

`git clone https://github.com/dngoldberg/MITgcm.git`

or

`git clone git@github.com:dngoldberg/MITgcm.git`

and checkout the branch:
`branch_runoff_2` (which is based on checkpoint 67v). Set your ROOTDIR environment variable to its location.

in `build/`, do the following:

`$ROOTDIR//tools/genmake2 -mods=../code -mpi -of=$ROOTDIR/tools/build_options/linux_amd64_gfortran --rootdir=$ROOTDIR`
`make depend`
`make -j 4`

### To prepare the initial ice-shelf state

open matlab, and run 

`runme`

### To run the experiment

Change line 5 of `runme.m` to `steps=4;`
In matlab, run `runme`

Note there are self-explanatory run time parameters at the top of the `'RunSingleCoupleStep'` section
