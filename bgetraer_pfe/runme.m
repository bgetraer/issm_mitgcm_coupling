% This runme controls a coupled MITgcm/ISSM model of a Pine Island Glacier-like ice shelf.
% The code has been adapted from https://github.com/MITgcm-contrib/issm_mitgcm_coupling.git
% 
% The runme has been altered to include the operations previously performed by external files such
% as input/rdmds_init.m, input/gendata.m, and ExpPar/SquareShelf.par. Additionally, the runme has
% been adapted to control the domain parameters in code/SIZE.h, input/data, and input/data.obcs.
% These changes allow for choices in domain resolution and size to be implemented consistently 
% from a single point of control.
% 
% Additionally, several changes have been made in how the model domain and boundary conditions are
% implemented, with the aim of more consistently reproducing the model described by Nakayama 2021
% and Jordan 2017.
%
% This file relies on matlab functions to read and write from SIZE.h, input/data, and input/data.obcs
% which are in ./m/
%
% The version here has been modified for running a coupled model on Pleiades.
% Some steps require MCC code to be compiled OUTSIDE of this runme script.
% The model run steps ('SteadystateNoSlip' and 'RunCouple') do not run within the MATLAB session,
% but run on the compute nodes by sending pre-compiled MCC executable MATLAB functions to the queue.


% Script control
steps=[3];		% organizer steps in this script to execute
addpath('./m');	% local matlab function directory
fbase=['/nobackup/bgetraer/issmjpl/proj-getraer/issm_mitgcm_coupling/bgetraer_pfe/'];	% project directory

devel=0; % send to devel queue or long queue?

% experiment control{{{
clear experiment;
experiment.name='QCTR'; % CTRL, INFL, GLCH, QCTR, QDST
% set CTRL options as default {{{
experiment.p=2; % inflow vel polynomial exponent
experiment.s=0; % geometry of grounding line induced channel (width=s, deltaH=0.15s) (m)
experiment.q=0; % total subglacial runoff flux (m^3/s)
experiment.runname=experiment.name;
% }}}
switch experiment.name
	case 'INFL' % Change pattern of inflow velocity {{{
		% p=(2),3,4,6
		experiment.p=6; % inflow vel polynomial exponent
		experiment.runname=sprintf('%s_p_%i',experiment.name,experiment.p);
		% }}}
	case 'GLCH' % Change geometry of grounding line induced channel {{{
		% s = 3000, 4333, 5667, 7000 m
		experiment.s=5667; % geometry of grounding line induced channel (width=s, deltaH=0.15s) (m)
		experiment.runname=sprintf('%s_s_%i',experiment.name,experiment.s);
		% }}}
	case 'QCTR' % Change subglacial runoff flux (point source in center) {{{
		% q=1,2,10,50 Gt/year
		experiment.q=0; % total subglacial runoff flux (Gt/year)
%		experiment.runname=sprintf('%s_q_%i',experiment.name,experiment.q);
%		experiment.s=4333; % geometry of grounding line induced channel (width=s, deltaH=0.15s) (m)
%     experiment.runname=sprintf('%s_s_%i',experiment.runname,experiment.s);
		experiment.w=3000; % width of grounding line induced channel (m)
      experiment.h=50;  % height of grounding line induced channel (m)
      experiment.runname=sprintf('%s_q_%i_w_%i_h_%i',experiment.name,experiment.q,experiment.w,experiment.h);
		% }}}
	case 'QDST' % Change subglacial runoff flux (distributed along inflow boundary) {{{
		% q=1,2,10,50 Gt/year
		experiment.q=50; % total subglacial runoff flux (Gt/year)
		experiment.runname=sprintf('%s_q_%i',experiment.name,experiment.q);
		% }}}
	case 'QCRF' % center flux runoff, refined ice near gr. line{{{
		experiment.w=3000; % width of grounding line induced channel (m)
		experiment.h=50;	% height of grounding line induced channel (m)
		experiment.q=0; % total subglacial runoff flux (Gt/year)
		experiment.runname=sprintf('%s_q_%i_w_%i_h_%i',experiment.name,experiment.q,experiment.w,experiment.h);
		% }}}
end
%set model name prefix 
prefix=sprintf('PigLike_%s_',experiment.runname);
% }}}
% define coordinate system {{{
	% define gridding coordinates {{{
	% The domain consists of a filled ocean domain surrounded on at least two sides by 
	% walls. The entire domain begins at (X0,Y0) with nWw cells between X0 and 0 and 
	% nWs cells between Y0 and 0. All cells outside of x=(0,LxOC) are walls, and all cells
	% outside of y=(0,LyOC) are walls.
	% The ISSM domain and MITgcm domains share a coordinate system, and cell centers and 
	% edges where they overlap.
	%
	% If changes are made to the coordinate system, 'BuildMITgcm' must be run to recompile
	% the model for the new domain.
	
	LxOC=60E3;  % length of ocean in x (m)
	LyOC=100E3; % length of ocean in y (m)
	Lz=1100;    % depth of ocean in z (m)
	LyICE=60E3;	% length of ice domain in y (m)
	dx=1e3;  % horizontal resolution in x (m)
	dy=dx;   % horizontal resolution in y (m)
	dz=20;   % vertical resolution in z (m)

	nWx=4; % number of wall cells in x
	nWy=2; % number of wall cells in y
	nWw=ceil(nWx/2);  % number of wall cells before x=0 (excess placed on east)
	nWs=nWy;          % number of wall cells before y=0 (excess placed on north)
	X0=0-nWw*dx; % origin of entire domain in x (m)
	Y0=0-nWs*dy; % origin of entire domain in y (m)
	
	NxOC=LxOC/dx;  % number of ocean cells in x
	NyOC=LyOC/dy;  % number of ocean cells in y
	Nx=NxOC+nWx;   % number of total cells in x
	Ny=NyOC+nWy;   % number of total cells in y
	Nz=Lz/dz;      % number of total cells in z

	xp=((0:Nx)*dx)+X0; % location of all cell edges in x (m)
	yp=((0:Ny)*dy)+Y0; % location of all cell edges in y (m)
	zp=(0:Nz)*dz;      % location of all cell edges in z (m)
   zg=zp(1:end-1); 	% location of z edge points, not including the lowest (redundant) boundary (m)
	
	xc=xp(1:end-1)+0.5*dx; % location of cell centers in x (m)
	yc=yp(1:end-1)+0.5*dy; % location of cell centers in y (m)
	zc=zp(1:end-1)+0.5*dz; % location of cell centers in z (m)

	xcOC=xc(xc>=0 & xc<=LxOC); % location of filled ocean cell centers in x (m)
	ycOC=yc(yc>=0 & yc<=LyOC); % location of filled ocean cell centers in y (m)

   % horizontal grid center points
   [XC YC]=meshgrid(xc,yc);
   XC=XC(:);
   YC=YC(:);
	indICE=(XC>=0 & XC<=LxOC & YC>=0 & YC<=LyICE); % XC and YC index for overlap with ice domain
	% }}}
	% define beta-plane for MITgcm coriolis {{{
	lat0=-75.5557; % reference latitude (deg)
	phi0=lat0/365*2*pi; % reference latitude (rad)
	rotationPeriod=8.6164E+04; % MITgcm default (s)
	radiusA=6.378E6; % equatorial radius (m)
	radiusB=6.357E6; % polar radius (m)
	radius0=1/sqrt((cos(phi0)/radiusA)^2 + (sin(phi0)/radiusB)^2);
	omega=2*pi/rotationPeriod; % angular velocity of Earth (rad/s)
	f0=2*omega*sin(phi0); % reference Coriolis parameter (1/s)
	beta0= 2*omega*cos(phi0)/radius0; % (1/(m*s))
	% write to MITgcm input/data {{{
	cd(fullfile(fbase,'input')); % move to input dir
	% write f0
	newline = [' f0 = ' num2str(round(f0,8)) ','];
   command=['sed "s/.*f0.*/' newline '/" data > data.temp; mv data.temp data'];
   system(command);
	% write beta
	newline = [' beta = ' num2str(round(beta0,14)) ','];
   command=['sed "s/.*beta.*/' newline '/" data > data.temp; mv data.temp data'];
   system(command);
	% selectCoriMap
	newline = [' selectCoriMap = 1,'];
   command=['sed "s/.*selectCoriMap.*/' newline '/" data > data.temp; mv data.temp data'];
   system(command);
	cd(fbase);
	% }}}
	% }}}
% }}}
% define coupled time-step parameters {{{
	coupled_time_step=1/365;	% length of coupled time step (y)
	nsteps=365*60;					% number of coupled time steps to take 
	MITgcmDeltaT=80;				% MITgcm time step (s)
	y2s=60*60*24*365;				% s/yr
	alpha_correction = 1;		% dmdt correction parameter: 1, fully "corrected", 0, no correction

	experiment.runname=sprintf('%s_dt_%i',experiment.runname,MITgcmDeltaT);
% }}}
% set cluster options {{{
cluster=generic('name',oshostname(),'np',27); % set number of processors for ISSM. 'name' will be filled at runtime
% }}}
% define experiment directory {{{
	expdir=fullfile(fbase,'experiments',experiment.runname);
	% make experiment directory if needed
	if ~exist(expdir)
		mkdir(expdir);
	end 
	% make model directory if needed
	modeldir=fullfile(expdir,'Models');
	if ~exist(modeldir)
      mkdir(modeldir);
   end
% }}}

org=organizer('repository',modeldir,'prefix',prefix,'steps',steps);
if perform(org,'BuildMITgcm'),% {{{
	% write domain parameters to input/data and input/data.obcs {{{
		writePARM04('./input/data',dx,dy,dz,Nx,Ny,Nz,'X0',X0,'Y0',Y0); % write PARM04 to input/data file
		writeDATAOBCS('./input/data.obcs',Nx,Ny); % write OB_Jnorth to input/data.obcs
	% }}}
	% define processing parameters and write to SIZE.h {{{
		% tiling must be adjusted for the number of total cells in the domain
		sNx=16;  % Number of X points in tile
		sNy=17;  % Number of Y points in tile
		OLx=3;   % Tile overlap extent in X
		OLy=3;   % Tile overlap extent in Y
		nSx=1;   % Number of tiles per process in X
		nSy=1;   % Number of tiles per process in Y
		nPx=Nx/sNx/nSx;   % Number of processes to use in X
		nPy=Ny/sNy/nSy;   % Number of processes to use in Y
		values=[sNx, sNy, OLx, OLy, nSx, nSy, nPx, nPy, Nx, Ny, Nz]; % default order (see writeSIZE.m)
		writeSIZE('./code/SIZE.h',values); % write to SIZE.h
	%}}}
   % recompile the MITgcm model {{{
      % MITgcm directories
      genmake2='${MITGCM_ROOTDIR}/tools/genmake2';
      optfile_dir='${MITGCM_ROOTDIR}/tools/build_options/linux_amd64_ifort+mpi_ice_nas';

      % clear the build directory
		cd(fbase);
		!rm -r ./build
		mkdir('build');
      cd('./build');
      % make the MITgcm executable
      command=[genmake2 ' -mpi -mo ../code -optfile ' optfile_dir ' -rd ${MITGCM_ROOTDIR}'];
  %    system(command); % generate Makefile
  % 	system('make CLEAN');	% prepare for new compilation
  %    system('make depend');  % create symbolic links from the local directory to the source file locations
  %    system('make');         % compile code and create executable file mitgcmuv
  %    cd(fbase);
   % }}}
end%}}}
if perform(org,'MeshParam'),% {{{
	% create mesh for ISSM model {{{
		md=model();	% initialize ISSM model structure
		%SQUAREMESH {{{
%		nnode_x=LxOC/dx+1;	% number of nodes in x
%		nnode_y=LyICE/dy+1;	% number of nodes in y
%		md=squaremesh(md,LxOC,LyICE,nnode_x,nnode_y); % make squaremesh
		% }}}
		%BAMG {{{
		% make exp file {{{
			corners=[0,(LxOC-dx)/2,LxOC,LxOC,0,0; 0,0,0,LyICE,LyICE,0];
			domainname='./Exp/domain.exp';
			fileID = fopen(domainname,'w');
			fprintf(fileID,'## Name:domainoutline\n');
			fprintf(fileID,'## Icon:0\n');
			fprintf(fileID,'# Points Count  Value\n');
			fprintf(fileID,'%i 1.0\n',size(corners,2));
			fprintf(fileID,'# X pos Y pos\n');
			fprintf(fileID,'%f %f\n',corners);
			fclose(fileID);
		% }}}
		%hmin=200; % minimum resolution along lateral boundaries (m)
		%hmax=750; % maximum resolution in interior of the model (m)
		%gradation=1.5;	% maximum ratio of element edge lengths
		md=triangle(model,domainname,dx); % initialize unstructured triangular mesh
		%hv=hmax*ones(md.mesh.numberofvertices,1); % initialize hvertices desired resolution (m)
		%%refine the corners of the ice shelf where we have large velocity gradient
		%ind=((md.mesh.x==0|md.mesh.x==LxOC)&(md.mesh.y==LyICE));
		%hv(ind)=hmin;	% set hvertices along the lateral edges (m)
		%md=bamg(md,'hVertices',hv,'gradation',gradation); % remesh
		if strcmp(experiment.name,'QCRF')
			hmin=200; % minimum resolution at center of ice (m)
			md=triangle(model,domainname,hmin); % initialize unstructured triangular mesh
			hmax=dx; % maximum resolution in interior of the model (m)
			gradation=1.5;   % maximum ratio of element edge lengths
			hv=hmax*ones(md.mesh.numberofvertices,1); % initialize hvertices desired resolution (m)
			hv(2)=hmin; % set hvertices at center inflow boundary (m)
			md=bamg(md,'hVertices',hv,'gradation',gradation); % remesh
		end
		% }}}
	% }}}
	% parameterize ISSM model and save {{{
	%  Set parameters and options for the ISSM model of the initial ice shelf
	%  before running uncoupled to steady state, and before channel initiation and so on.
	
	%Model name
	md.miscellaneous.name='PigLike';
	
	%Material propoerties
	md.materials.rheology_B=4.9E5*y2s^(1/3)*ones(md.mesh.numberofelements,1); % BENJY: changed to match Jordan, 2017 (Pa s^1/3)
	md.materials.rheology_n=3*ones(md.mesh.numberofelements,1); % Glen's Flow Law exponenent
	md.materials.rho_water=1028; % (kg/m^3)
	
	%Geometry
	%inflow thickness: if experiment.s>0 there is a channel with a gaussian shape at the center
	% at the grounding line this channel has a sigma of experiment.s/4,
	% i.e. experiment.s is equivalent to a width of 2sigma on either side.
	% the difference in ice thickness at the peak of the channel is given by lambda.
	H0=1200; % inflow thickness before perturbation (m)
	HLy=500; % thickness at edge of ice shelf
	if exist('experiment.h')
		lambda=experiment.h; % thickness difference at peak of grounding line induced channnel (m)
		experiment.s=experiment.w; % width of channel (m)
	else
		lambda=0.15*experiment.s; % thickness difference at peak of grounding line induced channnel (m)
	end
	x=md.mesh.x;y=md.mesh.y;
	Hx=H0-lambda.*exp(-1/2*(x-LxOC/2).^2 / (experiment.s/4).^2); % thickness in the x direction (m)
	if lambda==0; Hx=H0; end % account for singularity at experiment.s=0

	md.geometry.thickness=HLy*y/LyICE + Hx.*(1-y/LyICE); % initial thickness (m)
	md.geometry.base=-md.materials.rho_ice/md.materials.rho_water*md.geometry.thickness; % initial ice draft (m)
	md.geometry.surface=md.geometry.base+md.geometry.thickness; % initial ice surface (m)
	md.geometry.bed=-Lz*ones(md.mesh.numberofvertices,1); % bed depth (m)
	
	%Initial velocity and pressure
	md.initialization.vx=zeros(md.mesh.numberofvertices,1); % m/yr
	md.initialization.vy=zeros(md.mesh.numberofvertices,1); % m/yr
	md.initialization.vz=zeros(md.mesh.numberofvertices,1); % m/yr
	md.initialization.pressure=zeros(md.mesh.numberofvertices,1); % Pa
	
	%Friction 
	md.friction.coefficient=zeros(md.mesh.numberofvertices,1);
	md.friction.p=ones(md.mesh.numberofelements,1);
	md.friction.q=ones(md.mesh.numberofelements,1);
	
	%Temperature
	md.initialization.temperature=(273.15-15)*ones(md.mesh.numberofvertices,1); % deg K
	md.thermal.spctemperature=NaN*ones(md.mesh.numberofvertices,1);
	
	%Numerical parameters
	%Mass transport stabilization
	md.masstransport.stabilization=5; % Streamline upwind Petrov–Galerkin method (SUPG) 
	md.thermal.stabilization=1;
	md.stressbalance.restol=0.10;
	md.steadystate.reltol=0.02;
	md.stressbalance.reltol=0.02;
	md.stressbalance.abstol=NaN;
	md.transient.isthermal=0; %no thermodynamics within the ice
	md.transient.ismovingfront=0; %ice is lost at the constant front
	md.transient.isgroundingline=0; %no grounding line movement
	md.groundingline.migration='None';

	%Ice and ocean masks
	md.mask.ocean_levelset=-1*ones(md.mesh.numberofvertices,1); % ocean everywhere
	md.mask.ice_levelset=-1*ones(md.mesh.numberofvertices,1);	% ice everywhere
	md.mask.ice_levelset(md.mesh.y==LyICE)=1; % transition mask at edge of ice front
	
	% From SetIceShelfBC:
	%Free flow everywhere by default
	md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.referential=NaN*ones(md.mesh.numberofvertices,6); 
	md.stressbalance.loadingforce=0*ones(md.mesh.numberofvertices,3); 
	
	%Surface forcings
	md.smb.mass_balance=zeros(md.mesh.numberofvertices,1);
	
	%Basal forcings
	md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.geothermalflux=zeros(md.mesh.numberofvertices,1);
	
	%Thickness and Damage controls DO I NEED THIS?
	md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);
	md.balancethickness.thickening_rate=zeros(md.mesh.numberofvertices,1);
	md.balancethickness.spcthickness=NaN*ones(md.mesh.numberofvertices,1);
	md.damage.spcdamage=NaN*ones(md.mesh.numberofvertices,1);
	
	%Boundary conditions:
	%Note: BAMG mesh appears to place some nodes very close to, but not on the boundaries.
	%Ice front
	pos=find(md.mesh.y>LyICE-1); % ice front boundary
	md.stressbalance.spcvy(pos)=NaN; % no spcvy along the entire ice front
	%Lateral boundaries
	pos=find(md.mesh.x<1 | md.mesh.x>LxOC-1); % lateral boundaries
	md.stressbalance.spcvx(pos)=0; % no velocity into wall
	md.stressbalance.spcvy(pos)=0; % no slip along wall
	%Fixed flux inflow
	pos=find(md.mesh.y<1); % inflow boundary
	md.masstransport.spcthickness(pos)=md.geometry.thickness(pos); % keep the thickness constant
	b=LxOC/2;										% half width of the domain (m)
	%q=80E9;											% integrated flux across inflow boundary, per Nakayama, 2022 (m^3/yr)
	%Vc=q/(H0*4/3*b);								% velocity at center of domain, determined by flux (m/yr)
	Vc=2E3;											% velocity at center of domain, 2000 m/yr
	Vin=Vc*(1-abs(((md.mesh.x(pos)-b)./b).^experiment.p)); % the generalized equation for inflow velocity (m/yr)
	md.stressbalance.spcvy(pos)=Vin;			% fixed inflow velocity across boundary
	md.stressbalance.spcvx(pos)= 0;			% no velocity along boundary

	%Flow equation
	md=setflowequation(md,'SSA','all');

	%Save model
	savemodel(org,md);
	% }}}
md.cluster=cluster;
end%}}}
if perform(org,'SteadystateNoSlip'),% {{{
	md = loadmodel(org,'MeshParam');
	%set transient options
   md.timestepping = timesteppingadaptive(md.timestepping); % choose based on velocity
   md.timestepping.start_time=0;    % yr
   md.timestepping.final_time=200;  % yr
   md.settings.output_frequency=100;

   md.transient.requested_outputs={'default','BasalforcingsFloatingiceMeltingRate','Thickness','IceVolume'};
   md.verbose=verbose('convergence',0,'solution',1,'module',0);
   md.toolkits.DefaultAnalysis=bcgslbjacobioptions();

	md.cluster=cluster;

	%locations of org and queue script for runsteadystate
	orgfile=fullfile(org.repository,[org.prefix 'org.mat']);
	queuefile=fullfile(fbase,'runsteadystate','queue_runsteadystate.sh');
	%save the organizer and current model
	save(orgfile, 'org');
	savemodel(org,md);

	%execute queue_runsteadystate.sh with the current org
	system([queuefile ' ' org.repository ' ' orgfile]);
end%}}}
if perform(org,'RunCouple'),% {{{
	% build MITgcm experiment dir, set transient options {{{
		disp('************************************************************************************')
		cd(expdir); % move into experiment directory
		disp(['*   - building experiment dir ' pwd]);
		% replace expdir/input with bgetraer_pfe/input {{{
		if exist(fullfile(expdir,'input'))
			system(['rm -rf ' fullfile(expdir,'input')]);
		end
		system(['cp -r ' fullfile(fbase,'input') ' ' expdir]);
		% }}}
		cd(fullfile(expdir,'input')); % move into input directory
		disp(['*   - building input dir      ' pwd]);
		% build input directory {{{
		% define bathymetry; write to input file {{{
		   bathy=-Lz*ones(Ny,Nx); % m
		   walls=(XC<0 | XC>LxOC | YC<0); % index of wall locations
		   bathy(walls)=0; % m
		   fid = fopen('BATHY60.box','w','b'); fwrite(fid,bathy','real*8'); fclose(fid); % write to input file
		% }}}
		% define OBCS conditions; write to input files {{{
			% velocity profile at open boundary {{{
   		%  note: ensure that the volume flux is zero at the open boundary
   		   v_sfc=0.025;   % V at surface (m/s)
   		   v_bot=-0.025;  % V at bottom (m/s)
   		   v_ref=linspace(v_sfc,v_bot,Nz); % linear vertical velocity profile (m/s)
   		   V_OBC=repmat(v_ref,Nx,1); % matrix of velocities along the open boundary (m/s)
   		% }}}
   		% define thermo/halo cline parameters {{{
   		   z1=290; % depth of beginning of thermo/halo cline (m)
   		   z2=690; % depth of end of thermo/halo cline (m)
   		   cline=find(zc>=z1 & zc<=z2); % index of cells centers on the thermo/halo cline

   		   T_sfc=-1;   % potential temperature of the upper water mass (deg C)
   		   T_bot=1.2;  % potential temperature of the lower water mass (deg C)

   		   S_sfc=34;   % salinity of the upper water mass (g/kg)
   		   S_bot=34.7; % salinity of the lower water mass (g/kg)
   		% }}}
			% temperature profile at open boundary {{{
   		%  from Tmin at sfc increasing to Tmax at bottom
   		%  take Tmin and Tmax from the OBC of the original PIG experiment
   		   T_ref=nan(1,Nz);     % initialize vertical temp profile (deg C)
   		   T_ref(zc<z1)=T_sfc;  % upper water temp (deg C)
   		   T_ref(zc>z2)=T_bot;  % lower water temp (deg C)
   		   T_ref(cline)=linspace(T_sfc,T_bot,length(cline)); % linear temp profile along thermocline (deg C)
   		   T_OBC=repmat(T_ref,Nx,1); % matrix of temp along the open boundary (deg C)
   		% }}}
   		% salinity profile at open boundary {{{
   		%  from Smin at sfc increasing to Smax at bottom
   		%  take Smin and Smax from the OBC of the original PIG experiment
   		   S_ref=nan(1,Nz);     % initialize vertical salinity profile (g/kg)
   		   S_ref(zc<z1)=S_sfc;  % upper water salinity (g/kg)
   		   S_ref(zc>z2)=S_bot;  % lower water salinity (g/kg)
   		   S_ref(cline)=linspace(S_sfc,S_bot,length(cline)); % linear salinity profile along halocline (g/kg)
   		   S_OBC=repmat(S_ref,Nx,1); % matrix of salinity along the open boundary (g/kg)
   		% }}}
   		% write OBCS files {{{
   		   fid=fopen('vvel.obw','w','b');  fwrite(fid,V_OBC,'real*8'); fclose(fid);
   		   fid=fopen('theta.obw','w','b');  fwrite(fid,T_OBC,'real*8'); fclose(fid);
   		   fid=fopen('salt.obw','w','b');  fwrite(fid,S_OBC,'real*8'); fclose(fid);
   		% }}}
		% }}}
		% define initial T, S conditions; write to input files {{{
		   S_init=reshape(S_OBC,Nx,1,Nz); % take open boundary conditions for S (g/kg)
		   S_init=repmat(S_init,1,Ny,1);  % extend to all cells, no change in x (g/kg)
		   T_init=reshape(T_OBC,Nx,1,Nz); % take open boundary conditions for T (deg C)
		   T_init=repmat(T_init,1,Ny,1);  % extend to all cells, no change in x (deg C)
		
		   fid=fopen('theta.init','w','b');fwrite(fid,T_init,'real*8');fclose(fid);
		   fid=fopen('salt.init','w','b');fwrite(fid,S_init,'real*8');fclose(fid);
		% }}}
	% define initial ice shelf mass, open z cells, and free surface; write to input files {{{
	   % define material properties, constants, and equation of state {{{
	      rho_ice=917.;     % 917 is the default for SHELFICE and ISSM (kg/m^3)
	      g=9.81;           % m/s^2
	      rhoConst=1000.;   % density of freshwater (kg/m^3)
	      pa2db=1E-4;       % pascal to decibar
	      eos = 'jmd95z';   % equation of state from Jackett and McDougall (1995)
	   % }}}
		% read ice shelf thickness from ISSM {{{
	      md = loadmodel(org,'SteadystateNoSlip');
			VV=zeros(size(XC)); % initialize ice thickness from ISSM (m)
	      VV(indICE)=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(end).Thickness,XC(indICE),YC(indICE),'default',0); % m
	      VV = reshape(VV,[Ny Nx])'; % reshape from vector to matrix, take transpose bc MITgcm has i and j flipped. (m)
	      % BENJY: I dont think we use this file at all
	      binwrite('H_new.box',VV); % write to binary input file
	      shelficemass=VV*rho_ice; % ice mass per m^2 at cell centers (kg/m^2)
	   % }}}
		% calculate pressure, potential, etc. due to density {{{
	      p=zc*g*rhoConst*pa2db; % initial guess for pressure in decibar (1 db = 1E-4 kg/(s^2 m))
	      pcnvg=rms(p); % initialize convergence criterion
	      ptol=1e-13; % convergence tolerance for defining the pressure
	      %i=0; % benjy's loop counter, just for debugging
	      while pcnvg>ptol
	         p0=p;                                                    % save last pressure estimate (db)
	         rho=densjmd95(S_ref,T_ref,p);                            % get new estimate of density using equation of state (kg/m^3)
	         drho=rho-rhoConst;                                       % density anomaly (kg/m^3)
	         phiC=cumsum(dz*g*drho/rhoConst)-(dz/2)*g*drho/rhoConst;  % cumulative gravitational potential anomaly at the cell centers (m^2/s^2)
	         phiF=[0 cumsum(dz*g*drho/rhoConst)];                     % cumulative gravitational potential anomaly at ALL cell edges (m^2/s^2)
	         p=rhoConst*(zc*g+phiC)*pa2db;                            % new pressure estimate (db)
	         pcnvg=rms(p-p0);                                         % update convergence criterion
	         %i=i+1; disp(num2str(i)); % loop counter for debugging
	      end
	      massC=rhoConst*(phiC/g+ zc); % vertically integrated mass density at cell centers (kg/m^2)
	      massF=rhoConst*(phiF/g+ zp); % vertically integrated mass density at ALL cell edges (kg/m^2)
	   % }}}
		% calculate depth of the ice draft and free surface {{{
	      topo = zeros(Nx,Ny);    % initialize matrix for ice draft depth (m)
	      etainit = zeros(Nx,Ny); % initialize matrix for free surface
	      mergethreshold=0.29; % dimensionless merge threshold from Jordan, 2017
	      % loop through horizontal cells
	      for ix=1:Nx
	         for iy=1:Ny
	            mass=shelficemass(ix,iy);     % mass density of current cell (kg/m^2)
	            k = max(find(massF < mass));  % z cell where ice-ocean interface is (lowest cell where ice can displace water at upper cell edge)
	            if ~isempty(k) % there is at least SOME ice
	               if (mass < massC(k)) % ice CANNOT displace the water at the center of the kth cell
	                  topo(ix,iy)= -zg(k) - (mass-massF(k)) * (dz/2)/(massC(k)-massF(k));     % place the interface between the upper edge and the center
	               else                 % ice CAN displace the water at the center of the kth cell
	                  topo(ix,iy)= -zc(k) - (mass-massC(k)) * (dz/2)/(massF(k+1)-massC(k));   % place the interface between the center and the next edge
	               end
	               % round topo to the cell edge, set eta_init to make difference
	               dr=1 - (-zg(k)-topo(ix,iy))/dz; % fraction of kth cell that is open
	               if (dr > mergethreshold)
	                  % bring Ro_surf *up* to closest grid face & make etainit negative to compensate
	                  topo(ix,iy) = -zg(k);      % m
	                  etainit(ix,iy) = (dr-1)*dz;   % m
	               else
	                  % bring Ro_surf *down* to closest grid face & make etainit pos to compensate
	                  topo(ix,iy) = -zg(k+1); % m
	                  etainit(ix,iy) = (dr)*dz;     % m
	               end
	            end
	         end
	      end
	   % }}}
		% write to files {{{
	      fid = fopen('shelficemassinit.bin','w','b'); fwrite(fid,shelficemass,'real*8'); fclose(fid);
	      fid = fopen('shelftopo.round.bin','w','b'); fwrite(fid,topo,'real*8'); fclose(fid);
	      fid = fopen('etainit.round.bin','w','b'); fwrite(fid,etainit,'real*8'); fclose(fid);
	   % }}}
		% }}}
		% define runoff; write to input file {{{
		% BENJY: how does this work? what are the units?
		   flux = zeros(Nx,Ny);
			rho_fw=1E3; %density of freshwater kg/m^3
			switch experiment.name
				case {'QCTR', 'QCRF'}
					flux(NxOC/2+nWw,nWs+1) = experiment.q*1E12/rho_fw/y2s; % flux in m^3/s
				case 'QDST'
					flux(nWw+[1:NxOC],nWs+1) = experiment.q*1E12/rho_fw/y2s/NxOC; % flux in m^3/s
			end
			fid = fopen('runoff_flux.bin','w','b'); fwrite(fid,flux,'real*8'); fclose(fid);
		% }}}
		% }}}
      % rename previous run directory and create new one {{{
      if exist(fullfile(expdir,'run.old'))
          system(['\rm -rf ' fullfile(expdir,'run.old')]);
		end
      if exist(fullfile(expdir,'run'))
			system(['\mv ' fullfile(expdir,'run') ' ' fullfile(expdir,'run.old')]);
		end
		mkdir(fullfile(expdir,'run'));
		% }}}
		cd(fullfile(expdir,'run')); % move into run directory
		disp(['*   - building run dir        ' pwd]);
		% build run directory {{{
		!ln -s ../input/* .
		system(['ln -s ' fullfile(fbase,'build/mitgcmuv') ' .']);
		!rm eedata* data*
		!cp ../input/data* .
		!cp ../input/eedata .
		% }}}
		% set MIGgcm transient options {{{
      % set input/data options {{{
      newline = [' ntimesteps = ' num2str(coupled_time_step*y2s/MITgcmDeltaT + 1)];
      command=['!sed "s/.*ntimesteps.*/' newline '/" data > data.temp; mv data.temp data'];
      eval(command)
      newline = [' pChkptFreq = ' num2str(coupled_time_step*y2s)];
      command=['!sed "s/.*pChkptFreq.*/' newline '/" data > data.temp; mv data.temp data'];
      eval(command)
		newline = [' deltaT = ' num2str(MITgcmDeltaT)];
      command=['!sed "s/.*deltaT.*/' newline '/" data > data.temp; mv data.temp data'];
      eval(command)
		newline = [' monitorFreq=' num2str(MITgcmDeltaT)];
      command=['!sed "s/.*monitorFreq.*/' newline '/" data > data.temp; mv data.temp data'];
      eval(command)
      % }}}
      % set input/data.diagnostics options {{{
      newline = [' frequency(3) = ' num2str(coupled_time_step*y2s)];
      command=['!sed "s/.*frequency(3).*/' newline '/" data.diagnostics > data.temp; mv data.temp data.diagnostics'];
      eval(command)
      newline = [' frequency(4) = ' num2str(-coupled_time_step*y2s)];
      command=['!sed "s/.*frequency(4).*/' newline '/" data.diagnostics > data.temp; mv data.temp data.diagnostics'];
      eval(command)
      newline = [' timephase(4) = ' num2str(0)];
      command=['!sed "s/.*timephase(4).*/' newline '/" data.diagnostics > data.temp; mv data.temp data.diagnostics'];
      eval(command)
      % }}}
		% set input/data.shelfice options {{{
		switch experiment.name
			case {'QCTR', 'QDST','QCRF'}
				disp(['*   - SHELFICE runoff: yes']);
				newline=' SHELFICEaddrunoff=.TRUE.';
			otherwise
				disp(['*   - SHELFICE runoff: no']);
				newline=' SHELFICEaddrunoff=.FALSE.';
		end
		command=['!sed "s/.*SHELFICEaddrunoff.*/' newline '/" data.shelfice > data.temp; mv data.temp data.shelfice'];
		eval(command)
		% }}}
      n0=0;  % starting step (for the coupled loop)
		SZ=readSIZE([fbase 'code/SIZE.h']); % get the number of processors from SIZE.H	
		npMIT=SZ.values(7)*SZ.values(8); % number of processors for MITgcm
		% }}}
		cd(fbase); % move back to runme directory
		disp(['*   - returning to            ' pwd]);
	% }}}
	% load ISSM model and set transient options {{{
		disp(['*   - loading ISSM model    ' org.prefix 'SteadystateNoSlip']);
		% load and initialize ISSM model from results of uncoupled transient {{{
		md = loadmodel(org,'SteadystateNoSlip');
		% only keep the end of the steady state transient
		md.results.TransientSolution=md.results.TransientSolution(end); 
		md.geometry.base=md.results.TransientSolution.Base;
		md.geometry.thickness=md.results.TransientSolution.Thickness;
		md.geometry.surface=md.geometry.base+md.geometry.thickness;
		md.initialization.vx=md.results.TransientSolution.Vx;
		md.initialization.vy=md.results.TransientSolution.Vy;
		md.initialization.vel=md.results.TransientSolution.Vel;
		md.initialization.pressure=md.results.TransientSolution.Pressure;
		% }}}
		% set ISSM transient options {{{
		md.verbose=verbose('convergence',false,'solution',false,'control',false);
		md.timestepping=timestepping;
		md.timestepping.final_time=coupled_time_step;	% end of each ice run (yr)
		md.timestepping.time_step=1/365;	% length of ice time step (yr)
		md.timestepping.start_time=0;		% simulation starting time (yr)
		%md.transient.requested_outputs={'default','IceVolume'};
		md.settings.output_frequency=1;			

		md.miscellaneous.name=experiment.runname;
		% }}}
	% }}}
	% save env variables and execute runcouple on the queue {{{
		% declare the variables we want to save to envfile and pass to the MCC deployable
		vars={'Nx' 'Ny' 'NxOC' 'XC' 'YC' 'indICE'...                                    % coordinate system variables
			'coupled_time_step' 'n0' 'nsteps' 'MITgcmDeltaT' 'y2s' 'alpha_correction'... % timestepping variables
			'org' 'md' 'npMIT'}; % ISSM classes and number of processors for MITgcm
	   envfile=fullfile(org.repository,[org.prefix 'env.mat']); %locations env and queue script for runcouple
		
		disp(['*   - writing parameters to ' envfile]);
	   save(envfile,vars{:}); %save variables to envfile
	   %execute runcouple.m through mcc with the current env on the queue
		if devel
			queuefile=fullfile(fbase,'runcouple','develqueue_runcouple.sh');
		else
			%error('only devel rn');
			queuefile=fullfile(fbase,'runcouple','longqueue_runcouple.sh');
		end
		system([queuefile ' ' fullfile(expdir,'run') ' ' envfile]);
	% }}}
end%}}}

if perform(org,'PickupCouple'),% {{{
	disp('************************************************************************************')
	% timestepping parameters
	n0=10950; % the coupled step that we want to pickup from
	nsteps=365*70; % the step number we want to end at
	
	% load existing model
   prefix=org.prefix;	% temporarily change org prefix
   org.prefix=sprintf('%s%0.5i',prefix,n0);
	disp(['*   - loading model ' org.prefix 'RunCouple'])
   md=loadmodel(org,'RunCouple');
	%disp(['*   - loading model ' org.prefix 'PickupCouple'])
   %md=loadmodel(org,'PickupCouple');
   org.prefix=prefix;	% change prefix back

	% reinitialize model from pickup
	md.geometry.base=md.results.TransientSolution(end).Base;           % m
	md.geometry.thickness=md.results.TransientSolution(end).Thickness; % m
	md.geometry.surface=md.geometry.base+md.geometry.thickness;        % m
	md.initialization.vx=md.results.TransientSolution(end).Vx;         % m/yr
	md.initialization.vy=md.results.TransientSolution(end).Vy;         % m/yr
	md.initialization.vel=md.results.TransientSolution(end).Vel;       % m/yr
	md.initialization.pressure=md.results.TransientSolution(end).Pressure; % Pa

	% ensure name is unique
	md.miscellaneous.name=experiment.runname;
	
	% set npMIT
	SZ=readSIZE([fbase 'code/SIZE.h']); % get the number of processors from SIZE.H	
	npMIT=SZ.values(7)*SZ.values(8); % number of processors for MITgcm

	% declare the variables we want to save to envfile and pass to the MCC deployable
   vars={'Nx' 'Ny' 'NxOC' 'XC' 'YC' 'indICE'...                               % coordinate system variables
      'coupled_time_step' 'nsteps' 'MITgcmDeltaT' 'y2s' 'alpha_correction'... % timestepping variables
      'org' 'md' 'npMIT' 'n0'}; % ISSM classes and number of processors for MITgcm
   envfile=fullfile(org.repository,[org.prefix 'pickupenv.mat']); %locations env and queue script for runcouple
   disp(['*   - writing parameters to ' envfile]);
   save(envfile,vars{:}); %save variables to envfile

	%execute runcouple.m through mcc with the current env on the queue
      if devel
         queuefile=fullfile(fbase,'runcouple','develqueue_pickupcouple.sh');
      else
         %error('only devel rn');
         queuefile=fullfile(fbase,'runcouple','longqueue_pickupcouple.sh');
      end
      system([queuefile ' ' fullfile(expdir,'run') ' ' envfile]);
end%}}}


%['cp PigLike_CTRL_MeshParam.mat ' fullfile(org.repository,[org.prefix 'MeshParam.mat'])]
%['cp PigLike_CTRL_SteadystateNoSlip.mat ' fullfile(org.repository,[org.prefix 'SteadystateNoSlip.mat'])]
