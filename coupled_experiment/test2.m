alpha_correction=0;
nsteps=180;
steps=5;

eval('!cp Models/PigLike.RunCouple_0 Models/PigLike.RunCouple')
!rm run
!ln -s run_0 run

figure(1)
runme

figure(2)

alpha_correction=1;

eval('!cp Models/PigLike.RunCouple_1 Models/PigLike.RunCouple')
!rm run
!ln -s run_1 run
runme


