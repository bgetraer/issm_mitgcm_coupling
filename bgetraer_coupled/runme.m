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

% Script control
steps=[2:3];		% organizer steps in this script to execute
addpath('./m');	% local matlab function directory
fbase=[pwd '/'];	% project directory

% experiment control{{{
% CTRL (no induced channel)
% SG (initial conditions similar to Sergienko 2013, no-slip control experiment);
% GLCH or ROCH (grounding line or runoff induced channel)
experiment='CTRL_BAMG'; % CTRL (no induced channel), GLCH or ROCH (grounding line or runoff induced channel)
ch_gline=0;  % initiate grounding line channel
ch_runoff=0; % initiate runoff channel
% }}}
% define coordinate system {{{
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

	nWx=2; % number of wall cells in x
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
% define coupled time-step parameters {{{
	coupled_time_step=1/365;	% length of coupled time step (y)
	nsteps=365*5;					% number of coupled time steps to take 
	MITgcmDeltaT=100;				% MITgcm time step (s)
	y2s=60*60*24*365;				% s/yr
	alpha_correction = 0;		% dmdt correction parameter: 1, fully "corrected", 0, no correction
% }}}
% set cluster options {{{
clustername='totten';
switch clustername
   case 'totten'
      cluster=generic('name',oshostname(),'np',20);
   case 'pfe'
      cluster=pfe('cpuspernode',28,'numnodes',1,'time',20,'interactive',0,'processor','bro','queue','devel'); %time in minutes: devel/normal/long
   case 'amundsen'
      cluster=generic('name','amundsen.thayer.dartmouth.edu','np',20,'interactive',0);
   otherwise
      error('cluster not supported yet');
end
% }}}

org=organizer('repository',[fbase 'Models'],'prefix',['PigLike' experiment '_'],'steps',steps);

if perform(org,'BuildMITgcm'),% {{{
	% write domain parameters to input/data and input/data.obcs {{{
		writePARM04('./input/data',dx,dy,dz,Nx,Ny,Nz,'X0',X0,'Y0',Y0); % write PARM04 to input/data file
		writeDATAOBCS('./input/data.obcs',Nx,Ny); % write OB_Jnorth to input/data.obcs
	% }}}
	% define processing parameters and write to SIZE.h {{{
		% tiling must be adjusted for the number of total cells in the domain
		sNx=31;  % Number of X points in tile
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
		% MITgcm directories ($ROOTDIR points to /totten_1/bgetraer/MITgcm_dan/, Dan's fork of MITgcm)
		genmake2='$ROOTDIR/tools/genmake2';
		optfile_dir='$ROOTDIR/tools/build_options/linux_amd64_gfortran';
	
		% clear the build directory
		cd ./build
		!rm *
		% make the MITgcm executable
		%command=[genmake2 '-mpi -mo ../code -optfile ' optfile_dir ' -rd $ROOTDIR'];
		command=['!' genmake2 ' -mpi -mo ../code -rd $ROOTDIR'];
		eval(command); % generate Makefile
		eval('!make depend')    % create symbolic links from the local directory to the source file locations
		eval('!make')           % compile code and create executable file mitgcmuv
		cd ..
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
			corners=[0,LxOC,LxOC,0,0; 0,0,LyICE,LyICE,0];
			domainname='domain.exp';
			fileID = fopen(domainname,'w');
			fprintf(fileID,'## Name:domainoutline\n');
			fprintf(fileID,'## Icon:0\n');
			fprintf(fileID,'# Points Count  Value\n');
			fprintf(fileID,'%i 1.0\n',size(corners,2));
			fprintf(fileID,'# X pos Y pos\n');
			fprintf(fileID,'%f %f\n',corners);
			fclose(fileID);
		% }}}
		hmin=200; % minimum resolution along lateral boundaries (m)
		hmax=1000; % maximum resolution in interior of the model (m)
		gradation=1.5;	% maximum ratio of element edge lengths
		md=triangle(model,domainname,hmin); % initialize unstructured triangular mesh
		hv=hmax*ones(md.mesh.numberofvertices,1); % initialize hvertices desired resolution (m)
		%refine the corners of the ice shelf where we have large velocity gradient
		ind=((md.mesh.x==0|md.mesh.x==LxOC)&(md.mesh.y==LyICE));
		hv(ind)=hmin;	% set hvertices along the lateral edges (m)
		md=bamg(md,'hVertices',hv,'gradation',gradation); % remesh
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
	H0=1200; % thickness at the inflow boundary y=0 (m)
	%H=H0*ones(1,nnode_x); % thickness at inflow boundary (m)
	%if ch_gline
	%	cw=3E3; % channel width (m)
	%	ca=500E3; % channel area (m^2)
	%	H=H-normpdf(x,LxOC/2,cw)*ca; % thickness at inflow boundary (m)
	%end
	alpha=(500-H0)/LyICE; % slope of the initial ice shelf thickness equivalent to h(LyICE)=500m
	md.geometry.thickness=alpha*md.mesh.y+H0; % initial thickness (m)
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
	q=80E9;											% integrated flux across inflow boundary, per Nakayama, 2022 (m^3/yr)
	b=LxOC/2;										% half width of the domain (m)
	Vc=q/(H0*4/3*b);								% velocity at center of domain, determined by flux (m/yr)
	Vin=Vc*(1-((md.mesh.x(pos)-b)./b).^2);	% the parabolic equation for inflow velocity (m/yr)
	if strcmp(experiment,'SG')
		% properties for no-slip conditions from Sergienko, 2013
		Vc=1000;		% (m/yr)
		Vin=Vc*(1-((md.mesh.x(pos)-b)./b).^4);	% equation for inflow velocity (m/yr)
		Vin=Vc*(1-abs(((md.mesh.x(pos)-b)./b).^n)); % the parabolic equation for inflow velocity (m/yr)
	end
	md.stressbalance.spcvy(pos)=Vin;			% fixed inflow velocity across boundary
	md.stressbalance.spcvx(pos)= 0;			% no velocity along boundary

	%Flow equation
	md=setflowequation(md,'SSA','all');

	%Save model
	savemodel(org,md);
	% }}}
end%}}}
if perform(org,'SteadystateNoSlip'),% {{{
	md=loadmodel(org,'MeshParam');

	%transient model options
	md.timestepping = timesteppingadaptive(md.timestepping); % choose based on velocity
	md.timestepping.start_time=0;		% yr
	md.timestepping.final_time=200;	% yr
	md.settings.output_frequency=300; 

	md.transient.requested_outputs={'default','BasalforcingsFloatingiceMeltingRate','Thickness','IceVolume'};
	md.verbose=verbose('convergence',0,'solution',1,'module',0);
	md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
	md=solve(md,'Transient');

	savemodel(org,md);
end%}}}
if perform(org,'RunCouple'),% {{{
	% initialize MITgcm input files and build run directory {{{
		disp('************************************************************************************')
		disp('*	- initializing MITgcm files')
		% gendata {{{
		cd ./input/
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
		   flux(NxOC/2+nWw,nWs+1) = 0.0;
			fid = fopen('runoff_flux.bin','w','b'); fwrite(fid,flux,'real*8'); fclose(fid);
		% }}}
		cd ../
		% }}}
		% build run directory {{{
		disp('*	- building run directory')
      % rename previous run directory and create new one
      if exist ('run.old')
          !\rm -rf run.old
      end
      if exist ('run')
			!\mv run run.old
		end
		!\mkdir run
		cd ./run;
		!ln -s ../input/* .
		!ln -s ../build/mitgcmuv .
		!rm eedata* data*
		!cp ../input/data* .
		!cp ../input/eedata .
		% }}}
	% }}}
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
		%md.cluster.executionpath='/local/helene/issmjpl/proj-seroussi/TestCoupling/execution';
		%md.cluster.codepath = '/thayerfs/apps/issm/bin';
		%md.cluster.etcpath = '/thayerfs/apps/issm/etc';
		md.verbose=verbose('convergence',false,'solution',false,'control',false);
		md.timestepping.final_time=coupled_time_step;	% end of each ice run (yr)
		md.timestepping.time_step=1/365;	% length of ice time step (yr)
		md.timestepping.start_time=0;		% simulation starting time (yr)
		%md.transient.requested_outputs={'default','IceVolume'};
		md.settings.output_frequency=1;			
	% }}}
	% set MIGgcm transient options {{{
		t=0;	% current time (yr)
		newline = [' niter0 = ' num2str(t*y2s/MITgcmDeltaT)];
   	command=['!sed "s/.*niter0.*/' newline '/" data > data.temp; mv data.temp data'];
   	eval(command)
   	newline = [' ntimesteps = ' num2str(coupled_time_step*y2s/MITgcmDeltaT + 1)];
   	command=['!sed "s/.*ntimesteps.*/' newline '/" data > data.temp; mv data.temp data'];
   	eval(command)
   	newline = [' frequency(3) = ' num2str(coupled_time_step*y2s)];
   	command=['!sed "s/.*frequency(3).*/' newline '/" data.diagnostics > data.temp; mv data.temp data.diagnostics'];
   	eval(command)
   	newline = [' frequency(4) = ' num2str(-coupled_time_step*y2s)];
   	command=['!sed "s/.*frequency(4).*/' newline '/" data.diagnostics > data.temp; mv data.temp data.diagnostics'];
   	eval(command)
   	newline = [' timephase(4) = ' num2str(0)];
   	command=['!sed "s/.*timephase(4).*/' newline '/" data.diagnostics > data.temp; mv data.temp data.diagnostics'];
   	eval(command)
   	newline = [' pChkptFreq = ' num2str(coupled_time_step*y2s)];
   	command=['!sed "s/.*pChkptFreq.*/' newline '/" data > data.temp; mv data.temp data'];
   	eval(command)
	% }}}
	% loop through each coupled step, run the models, save the ouput {{{
		disp('*  - begin coupled model')
		SZ=readSIZE([fbase 'code/SIZE.h']); % get the number of processors from SIZE.H	
		npMIT=SZ.values(7)*SZ.values(8); % number of processors for MITgcm
		results=md.results;	% initialize results structure to store ISSM results (keep the end of the steadystate transient)
		% BEGIN THE LOOP
		%	n is the number step we are on, from 0:nsteps-1
		for n=0:(nsteps-1);
			display(['STEP ' num2str(n) '/' num2str(nsteps-1)]);
			t=n*coupled_time_step; % current time (yr)
			% calculate difference between ISSM and MITgcm mass, and corresponding adjustment to dmdt {{{
			if (n>0);
				mitgcm_mass=rdmds([fbase 'run/SHICE_mass'], round(t*y2s/MITgcmDeltaT)); % MITgcm mass at cell centers (kg/m^2)
				VV=zeros(size(XC)); % initialize ice thickness from ISSM (m)
		      VV(indICE)=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(end).Thickness,XC(indICE),YC(indICE),'default',0); % m
		      VV = reshape(VV,[Ny Nx])'; % reshape from vector to matrix, take transpose bc MITgcm has i and j flipped. (m)
				issm_mass=VV*rho_ice; % ice mass per m^2 at cell centers (kg/m^2)
				dmdt_adjust=alpha_correction*(issm_mass-mitgcm_mass)/(coupled_time_step*y2s); % difference per coupled time step (kg/m^2/s) 
				%dmdt_adjust(mitgcm_mass==0)=0;
			else
				dmdt_adjust=zeros(Ny,Nx)';
			end
			% }}}
			% RUN ISSM without melt{{{
				md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1); % ensure zero melt
				tic
				md=solve(md,'Transient'); 
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

				binwrite([fbase 'run/shelfice_dmdt.bin'],dmdt + dmdt_adjust);
				system(['cp ' fbase 'run/shelfice_dmdt.bin ' fbase 'run/shelfice_dmdt_' num2str(n) '.bin']);
			% }}}
			% RUN MITgcm {{{
				% update MITgcm transient options 
				newline = [' niter0 = ' num2str(t*y2s/MITgcmDeltaT)];
				command=['!sed "s/.*niter0.*/' newline '/" data > data.temp; mv data.temp data'];
				eval(command)
				% run MITgcm
				disp('about to run MITgcm')
				tic
				eval(['!mpirun -np ' int2str(npMIT) ' ./mitgcmuv > out 2> err']);
				disp('done MITgcm')
				toc
			% }}}
			% read melt from MITgcm run and write to ISSM model {{{ 
				system(['cp ' fbase 'run/SHICE_fwFluxtave.' appNum(round((t+coupled_time_step)*y2s/MITgcmDeltaT),10) '.data ' fbase 'run/melt.data']);
				system(['cp ' fbase 'run/STDOUT.0000 ' fbase 'run/stdout' num2str(n)]);

				melt_mitgcm=binread('melt.data',4,Nx,Ny)'; % melt flux at cell centers (kg/m^2/s)
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
				md=solve(md,'Transient');
			% }}}
			% sav results of ISSM run and reinitialize model {{{
				results.TransientSolution(end+1)= md.results.TransientSolution(end);
				results.TransientSolution(end).time = (n+1)*coupled_time_step;

				md.geometry.base=md.results.TransientSolution(end).Base;				 % m
				md.geometry.thickness=md.results.TransientSolution(end).Thickness; % m
				md.geometry.surface=md.geometry.base+md.geometry.thickness;			 % m
				md.initialization.vx=md.results.TransientSolution(end).Vx;			 % m/yr
				md.initialization.vy=md.results.TransientSolution(end).Vy;			 % m/yr
				md.initialization.vel=md.results.TransientSolution(end).Vel;		 % m/yr
				md.initialization.pressure=md.results.TransientSolution(end).Pressure; % Pa
			% }}}
		end
	% }}}
	% save ISSM model {{{
		md.results = results;
		savemodel(org,md);
		cd ..
	% }}}
end%}}}
if perform(org,'PlotOutput') % {{{ 

   md = loadmodel(org,'RunCouple');

	filepath = './run/';
	%calculate the locations for the cell corners from the cell centers
	xg=xc-0.5*dx;
	yg=YC-0.5*dy;
	% U points (on the Arakawa C grid) are on XG, YC
	% V points (on the Arakawa C grid) are on XC, YG


	t_file=round((nsteps-1)*y2s*coupled_time_step/MITgcmDeltaT);
	shelficethick = rdmds('run/SHICE_mass',t_file)/md.materials.rho_ice;
	ocean_v = rdmds('run/V',t_file);
	ocean_u = rdmds('run/U',t_file);
	%read in the free surface height after 3 years (77760 steps * 1200 seconds) from the ouput files
	Eta=rdmds([filepath 'Eta'], t_file);
	d = rdmds([filepath 'Depth']);	
	%plot
	figure(1);clf;hold on;
	subplot(3,1,1)
	imagesc(yg,xc,ocean_v(:,:,30));
	axis equal tight
	set(gca,'xticklabel',get(gca,'xtick')*1E-3)
	set(gca,'yticklabel',get(gca,'ytick')*1E-3)
	xlabel('km'); ylabel('km');
	title('MITgcm ocean bottom V');
	cb = colorbar;
   cb.Label.String = 'Velocity (m/s)';

	subplot(3,1,2)
   imagesc(yc,xg,ocean_u(:,:,30));
   axis equal tight
	set(gca,'xticklabel',get(gca,'xtick')*1E-3)
	set(gca,'yticklabel',get(gca,'ytick')*1E-3)
   xlabel('km'); ylabel('km');
   title('MITgcm ocean bottom U');
   cb = colorbar;
   cb.Label.String = 'Velocity (m/s)';
	
	subplot(3,1,3)
	A = shelficethick-mean(shelficethick(2:end,:));
	A(shelficethick<1)=0;
	imagesc(yc,xc,shelficethick)
	axis equal tight
	set(gca,'xticklabel',get(gca,'xtick')*1E-3)
	set(gca,'yticklabel',get(gca,'ytick')*1E-3)
	xlabel('km'); ylabel('km');
	title('MITgcm ice draft');
	cb = colorbar;
	cb.Label.String = 'Depth (m)';

end % }}}

if perform(org,'TestDrift'),% {{{ 
	md = loadmodel(org,'RunCouple');
	for i=0:(nsteps-1);

		if (i>0);
			shelficethick = rdmds('run/SHICE_mass',round(i*y2s*coupled_time_step/MITgcmDeltaT))'/md.materials.rho_ice;
		else
			shelficethick = binread('run/shelficemassinit.bin',8,60,100)'/md.materials.rho_ice;
		end
		if (i==0);
			shelficethick0=shelficethick;
		end
		issm_thick = md.results.TransientSolution(i+1).Thickness;
		shelficethickmesh = InterpFromGridToMesh(xc',yc',shelficethick,md.mesh.x,md.mesh.y,0);
		issm_thick_mitgcm = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,issm_thick,XC,YC,'default',0);
		issm_thick_mitgcm = reshape(issm_thick_mitgcm,[Ny Nx]);
		if (i==0);
			issm_thick_mitgcm0 = issm_thick_mitgcm;
		end

		diff = shelficethickmesh-issm_thick;
		diff(md.mesh.x<1.1e3 | md.mesh.y<1.1e3) = 0;
		diff_grid = shelficethick - issm_thick_mitgcm;
		diff_grid(:,1) = 0;
		diff_grid(1,:) = 0;
	end

	subplot(1,3,1);
	rge = prctile(abs(issm_thick_mitgcm(:)-issm_thick_mitgcm0(:)),99);
	pcolor(issm_thick_mitgcm-issm_thick_mitgcm0); colorbar; shading flat; caxis([-rge rge]);
	subplot(1,3,2); 
	pcolor(shelficethick-shelficethick0); colorbar; shading flat; caxis([-rge rge]);
	subplot(1,3,3); 
	diff=(shelficethick-shelficethick0)-(issm_thick_mitgcm-issm_thick_mitgcm0);
	diff(1,:)=0;
	diff(:,[1 end])=0;
	pcolor(diff); colorbar; shading flat;

end%}}}
