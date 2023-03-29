devpath
alpha_correction=0;
nsteps=180;
steps=4;
runme

eval('!cp -r run run_0')
eval('!cp Models/PigLike.RunCouple Models/PigLike.RunCouple_0')



alpha_correction=1;
nsteps=180;
steps=4;
runme

eval('!cp -r run run_1')
eval('!cp Models/PigLike.RunCouple Models/PigLike.RunCouple_1')
