function runsteadystate(orgfile)
%RUNSTEADYSTATE is intended as a barebones script to run a transient ISSM solve with MCC compilation
%The only input is the "org" structure, which must contain the step 'SteadystateNoSlip.' 
%RUNSTEADYSTATE loads the existing model with transient model options, solves, and saves.
%
%DEPLOYMENT:
%	IN MATLAB
%  createMCC('runsteadystate.m')
%	IN BASH
%  ./mccfiles/run_MCCexecutable.sh $ISSM_DIR/lib:$PETSC_DIR/lib:$MPI_ROOT/lib:${MKLROOT}/lib/intel64_lin:${MKLROOT}/../compiler/lib/intel64_lin:${ISSM_DIR}/externalpackages/triangle/install/lib:/nasa/matlab/2022b orgfile

%load the organizer
org=organizer();
load(orgfile);

%load SteadystateNoSlip model from runme
md=model();
md.timestepping=timesteppingadaptive();
md=loadmodel(org,'SteadystateNoSlip');

%solve and save
md.cluster.name=oshostname();
md=solve(md,'Transient');
savemodel(org,md);
