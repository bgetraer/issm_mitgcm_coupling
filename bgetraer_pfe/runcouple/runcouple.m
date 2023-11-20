function runcouple(envfile)
%RUNCOUPLE is intended as a standalone script to run a coupled ISSM-MITGCM model with MCC compilation/
%The only input is envfile.mat, which must contain all environment variables needed for this function,
%exported from the main runme.m. These variables are declared explicitly and named in this script as 
%the executable needs to expect them in order to load them.
%RUNCOUPLE loads the existing environment variables, loops through the time steps calling the models runs,
%and saves the output. Is is assumed that you are already located within the mitgcm "run" directory. 

dispMITxISSM();
disp('************************************************************************************');
disp('*   - beginning RUNCOUPLE mcc deployable');
disp(['*   - current directory is ' pwd])
%declare variables and load environment {{{
	%declare all variables we need to load from envfile {{{
	% declare coordinate system variables
	vars={'Nx' 'Ny' 'NxOC' 'XC' 'YC' 'indICE'};
	Nx=[];
	Ny=[];
	NxOC=[];
	XC=[];
	YC=[];
	indICE=[];
	% declare timestepping variables
	vars={vars{:} 'coupled_time_step' 'nsteps' 'MITgcmDeltaT' 'y2s' 'alpha_correction' 'n0'};
	coupled_time_step=[];
   nsteps=[];
   MITgcmDeltaT=[];
   y2s=[];
   alpha_correction=[];
	n0=[];
	% declare ISSM specific classes 
	vars={vars{:} 'org' 'md'};
   org=organizer();
	md=model();
	% declare MITgcm processors
	vars={vars{:} 'npMIT'};
	npMIT=[];
	% }}}
	%load environment variables and initialize ISSM cluster and results {{{
	load(envfile,vars{:});         % load environment variables
	md.cluster.name=oshostname();  % set cluster name
	results=md.results;            % initialize results structure to store ISSM results (keep the end of the steadystate transient)
	% }}}
% }}}
%loop through each coupled step, run the models, save the ouput {{{
	% n is the number step we are on, from 0:nsteps-1
	for n=n0:(nsteps-1);
		display(['STEP ' num2str(n) '/' num2str(nsteps-1)]);
		t=n*coupled_time_step; % current time (yr)
		% set when to write output to permanent files {{{
			% save md every 6 months and at the end
			writeISSM=(any(mod(n+1,365)==[0,180]) | n==(nsteps-1));
			% keep pickup files once every 30 days and at the end of the year
			% last step is kept until next step finishes.
         saveMITgcm=any(mod(n,365)==[0:30:330]);
		% }}}
		% calculate difference between ISSM and MITgcm mass, and corresponding adjustment to dmdt {{{
		if (n>0);
			filename=sprintf('SHICE_mass.%0.10i.data',round(t*y2s/MITgcmDeltaT));
			mitgcm_mass=binread(filename,4,[Nx,Ny]); % MITgcm mass at cell centers (kg/m^2)
			VV=zeros(size(XC)); % initialize ice thickness from ISSM (m)
	      VV(indICE)=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(end).Thickness,XC(indICE),YC(indICE),'default',0); % m
	      VV = reshape(VV,[Ny Nx])'; % reshape from vector to matrix, take transpose bc MITgcm has i and j flipped. (m)
			issm_mass=VV*md.materials.rho_ice; % ice mass per m^2 at cell centers (kg/m^2)
			dmdt_adjust=alpha_correction*(issm_mass-mitgcm_mass)/(coupled_time_step*y2s); % difference per coupled time step (kg/m^2/s) 
			%dmdt_adjust(mitgcm_mass==0)=0;
		else
			dmdt_adjust=zeros(Ny,Nx)';
		end
		% }}}
		% RUN ISSM without melt{{{
			md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1); % ensure zero melt
			tic
			md=solve(md,'Transient','runtimename',0);
			toc
		% }}}
		% calculate dmdt from ISSM run and write to MITgcm file {{{
			dH_issm=(md.results.TransientSolution(end).Thickness-md.geometry.thickness);  % change in thickness at element vertices (m)
			dt_issm=(md.timestepping.final_time-md.timestepping.start_time);					% ISSM time interval (yr)
			dmdt_issm=md.materials.rho_ice*dH_issm/dt_issm;											% dmdt at element vertices (kg/m^2/yr)
			dmdt_issm(md.results.TransientSolution(end).Thickness<=2)=0;						% dmdt=0 where there is "no ice" (kg/m^2/yr)
			% interpolate dmdt from ice coordinates to cell centers
			dmdt=zeros(size(XC)); % initialize dmdt at cell centers (kg/m^2/yr)
			dmdt(indICE)=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,dmdt_issm,XC(indICE),YC(indICE),'default',0); % (kg/m^2/yr)
			dmdt=1/y2s*reshape(dmdt,[Ny,Nx])'; % dmdt at cell centers, reshaped and in units for MITgcm (kg/m^2/s)

			% write to MITgcm read file (overwritten every step)
			binwrite('shelfice_dmdt.bin',dmdt+dmdt_adjust,8);
			% niter of the end of this run (will save every saveMITgcm)
			thisniter=round(t+coupled_time_step*y2s/MITgcmDeltaT);
			% write to numbered file to track
			binwrite(sprintf('shelfice_dmdt_%0.10i.bin',thisniter),dmdt+dmdt_adjust,8);
		% }}}
		% RUN MITgcm {{{
			% update MITgcm transient options
			newline = [' niter0 = ' num2str(t*y2s/MITgcmDeltaT)];
			command=['sed "s/.*niter0.*/' newline '/" data > data.temp; mv data.temp data'];
			system(command);
			
			% run MITgcm
			disp('about to run MITgcm')
			tic
			system(['mpirun -np ' int2str(npMIT) ' ./mitgcmuv > out 2> err']);
			disp('done MITgcm')
			toc

			% if model crashed, save ISSM state and end script
			errlines=readlines('err'); % read err file
			if any(contains(errlines,'ABNORMAL END'))
				disp(sprintf('ABNORMAL END to MITgcm on coupling step %i/%i',n,nsteps-1));
				disp('Saving state of ISSM model:');
				% save ISSM model
				md.results.TransientSolution(1).time=(n+1)*coupled_time_step;
				org.prefix=sprintf('%s%0.5i',org.prefix,n+1);
				savemodel(org,md);
				% throw error
				error('Ending RunCouple due to model crash');
			end
		% }}}
		% read melt from MITgcm run and write to ISSM model {{{ 
			% niter of this finished run
			thisniter=round((t+coupled_time_step)*y2s/MITgcmDeltaT);
			% save MITgcm output to numbered file
			filename=sprintf('stdout.%0.10i',thisniter);
			system(['cp ./STDOUT.0000 ' filename]);

			% locate the melt output of MITgcm
			filename=sprintf('SHICE_fwFluxtave.%0.10i.data',thisniter);
			melt_mitgcm=binread(filename,4,[Nx,Ny])'; % melt flux at cell centers (kg/m^2/s)
			melt_mitgcm=reshape(melt_mitgcm,[Ny,Nx]);  % (kg/m^2/s)
			% interpolate/extrpolate melt from MITgcm cell centers onto ISSM nodes
			X1=reshape(XC(indICE),NxOC,NxOC); % x location of cell centers within ice domain (m)
			X2=reshape(YC(indICE),NxOC,NxOC); % y location of cell centers within ice domain (m)
			V=reshape(melt_mitgcm(indICE),NxOC,NxOC); % melt flux only within the ice covered domain
			F=griddedInterpolant(X1',X2',V','linear','linear'); % linear interpolant and linear extrapolant
			%Ice mesh center of elements
			xe = mean(md.mesh.x(md.mesh.elements),2); ye = mean(md.mesh.y(md.mesh.elements),2);
			%meltq=F(md.mesh.x,md.mesh.y); % the melt at the queried ISSM nodes (kg/m^2/s)
			meltq=F(xe,ye); % the melt at the queried ISSM nodes (kg/m^2/s)
			%meltq(md.mesh.x==0 | md.mesh.x==LxOC | md.mesh.y==0)=0; % enforce zero melt along the no-slip and fixed inflow boundaries

			md.basalforcings.floatingice_melting_rate=-meltq*y2s/md.materials.rho_ice; % melt rate at element vertices (m/yr)
		% }}}
		% RUN ISSM with melt {{{
			md=solve(md,'Transient','runtimename',0);
		% }}}
		% save results of ISSM run and reinitialize model {{{
			%save md every 6 months and at the end
			if (writeISSM)
				% store most recent TransientSolution in temp structure
				mdtemp=md;
				mdtemp.results.TransientSolution=md.results.TransientSolution(end);
				mdtemp.results.TransientSolution(1).time=(n+1)*coupled_time_step;
				% temporarily change org prefix to save model
				prefix=org.prefix;
				org.prefix=sprintf('%s%0.5i',prefix,n+1);
				savemodel(org,mdtemp);
				org.prefix=prefix;
			end
			
			%reinitialize model
			md.geometry.base=md.results.TransientSolution(end).Base;				 % m
			md.geometry.thickness=md.results.TransientSolution(end).Thickness; % m
			md.geometry.surface=md.geometry.base+md.geometry.thickness;			 % m
			md.initialization.vx=md.results.TransientSolution(end).Vx;			 % m/yr
			md.initialization.vy=md.results.TransientSolution(end).Vy;			 % m/yr
			md.initialization.vel=md.results.TransientSolution(end).Vel;		 % m/yr
			md.initialization.pressure=md.results.TransientSolution(end).Pressure; % Pa
		% }}}
		 % remove unneeded MITgcm files {{{
         % n is the current step, representing the END OF THE PREVIOUS time step
         prev_end_niter=round(t*y2s/MITgcmDeltaT); % end of the PREVIOUS time step
         this_beg_niter=round(t*y2s/MITgcmDeltaT+1); % beginning of THIS time step
         if saveMITgcm & n>n0
            % do not delete end of last time step
            command=sprintf('rm *%0.10i*',this_beg_niter);
         elseif n>n0
            % delete end of last time step
            command=sprintf('rm *%0.10i* *%0.10i*',this_beg_niter,prev_end_niter);
         end
         system(command);
         % }}}
	end
% }}}
disp('************************************************************************************');
disp('*   - RUNCOUPLE finished');
disp('************************************************************************************');

% subfunctions
function D=binread(fname,prec,arrsize) % {{{
% read data from binary file into a matlab array D.
% Assumes big-endian architecture, and given precision
% and array size.
%
% fname: filename or path (string)
% prec: 4 or 8 for number of bits
% arrsize: dimensions of D (array)
%
% D: array of requested dimension
	fid=fopen(fname,'r','b');
	switch prec
		case 8
			D=fread(fid,inf,'real*8');
		case 4
			D=fread(fid,inf,'real*4');
		otherwise
			error('give precision of data');
	end
	D=reshape(D,arrsize);
	fclose(fid);
end % }}}
function q=binwrite(fname,D,prec) % {{{
% write a matlab array D of arbitrary dimension to binary file
% using big-endian architecture and given precision.
%
% fname: filename or path (string)
% D: array of arbitrary dimension (storage is independent of dimension
% sizes)
% prec: 4 or 8 for number of bits
	fid=fopen(fname,'w','b');
	switch prec
		case 8
			q=fwrite(fid,D,'real*8');
		case 4
			q=fwrite(fid,D,'real*4');
		otherwise
			error('use valid precision');
	end
	fclose(fid);
end % }}}
function dispMITxISSM() % {{{
%DISPMITXISSM prints an ascii graphic to terminal output
   disp(' __      __   _   _______                    _     _  _____   _____ ____ __      __ ')
   disp('|   \  /   | | | |__   __|                   \ \ / / |_   _|/ ____/ ____|   \  /   |')
   disp('| |\ \/ /| | | |    | | __ _  ___ _ __ ___    \   /    | |  \___ \\___ \| |\ \/ /| |')
   disp('| | \__/ | | | |    | |/ _` |/ __| `_ ` _ \   /   \   _| |_ ____| |___| | | \__/ | |')
   disp('|_|      |_| |_|    |_| (_| | (__| | | | | | /_/ \_\ |_____|_____/_____/|_|      |_|')
   disp('                       \___ |\___|_| |_| |_| ')
   disp('                        __/ |                ')
   disp('                       |___/                 ')
end % }}}
end
