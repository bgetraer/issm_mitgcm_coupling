function gatherResults(expname,startyear,endyear)
%GATHERRESULTS reads results from ISSM and MITGCM output and saves them each to a .mat file
% Useage:
%   gatherResults(expname,startyear,endyear)
%      where expname is of the form 'CTRL_dt_100' for example.
%      results include everthing AFTER the startyear, until and including the endyear.
%      i.e. gatherResults(expname,1,4) will return results for days 366 through 365*4 

% define repositories
exprepo='/nobackupp11/bgetraer/issmjpl/proj-getraer/issm_mitgcm_coupling/bgetraer_pfe/experiments';
expdir=fullfile(exprepo,expname); % where these results are located
resultsrepo='/nobackupp11/bgetraer/issmjpl/proj-getraer/issm_mitgcm_coupling/bgetraer_pfe/results';

% file indexing by coupled run step {{{
%startyear=0; % start time in years
%endyear=40; % end time in years
alln=(startyear*365):(endyear*365); % all n that we looped through (coupled step or day)
alln=alln(alln>startyear*365); % ignore day 0 of the start year

nstepissm=find(any(mod(alln,365)==[0,180]')); % steps where we saved ISSM
nstepmitgcm=find(any(mod(alln,365)==[0:30:330]')); % steps where we saved MITgcm
nissm=length(nstepissm);
nmitgcm=length(nstepmitgcm);
% }}}
% gather ISSM results {{{
% issm coordinates {{{
dx=1000; %m
Lx=60E3; %m
Ly=60E3; %m
xc=dx/2:dx:Lx;%cell center in x (m)
yc=dx/2:dx:Ly;%cell center in y (m)
[XCice,YCice]=meshgrid(xc,yc);
% }}}
% initialize ISSM data structures {{{
nvertices=2938; % number of ISSM vertices
nelements=5668; % number of ISSM elements

issmVx=zeros(nvertices,nissm);% ice Vx m/yr
issmVy=zeros(nvertices,nissm);% ice Vy m/yr
issmVel=zeros(nvertices,nissm);% ice Vel m/yr
issmThickness=zeros(nvertices,nissm);% ice thickness m
issmBasalmeltrate=zeros(nelements,nissm);% basal melt rate m/yr

interpVx=zeros(length(yc),length(xc),nissm);
interpVy=zeros(length(yc),length(xc),nissm);
interpVel=zeros(length(yc),length(xc),nissm);
interpThickness=zeros(length(yc),length(xc),nissm);
interpBasalmeltrate=zeros(length(yc),length(xc),nissm);
% }}}
disp(['Loading ISSM results from     ' fullfile(expdir,'Models')]);
% gather ISSM results {{{
for i=1:nissm
	disp(sprintf('reading step %i/%i',i,nissm));
	%load model
	filename=sprintf('*%0.5i*Couple.mat',nstepissm(i));
	filepath=dir(fullfile(expdir,'Models',filename));
	md=loadmodel(fullfile(filepath.folder,filepath.name));
	%save results on ISSM vertices and elements
	issmVx(:,i)=md.results.TransientSolution(end).Vx; % ice Vx m/yr
	issmVy(:,i)=md.results.TransientSolution(end).Vy; % ice Vy m/yr
	issmVel(:,i)=md.results.TransientSolution(end).Vel; % ice Vel m/yr
	issmThickness(:,i)=md.results.TransientSolution(end).Thickness; % ice thickness m
	issmBasalmeltrate(:,i)=md.results.TransientSolution(end).BasalforcingsFloatingiceMeltingRate; % basal melt rate m/yr	
	% interpolated onto MITgcm cell center grid
	interpVx(:,:,i)=griddata(md.mesh.x,md.mesh.y,issmVx(:,i),XCice,YCice);
	interpVy(:,:,i)=griddata(md.mesh.x,md.mesh.y,issmVy(:,i),XCice,YCice);
	interpVel(:,:,i)=griddata(md.mesh.x,md.mesh.y,issmVel(:,i),XCice,YCice);
	interpThickness(:,:,i)=griddata(md.mesh.x,md.mesh.y,issmThickness(:,i),XCice,YCice);
	elX=mean(md.mesh.x(md.mesh.elements),2); elY=mean(md.mesh.y(md.mesh.elements),2);
	interpBasalmeltrate(:,:,i)=griddata(elX,elY,issmBasalmeltrate(:,i),XCice,YCice);
end
issmX=md.mesh.x; % ice x coordinates (m)
issmY=md.mesh.y; % ice y coordinates (m)
interpX=XCice;
interpY=YCice;
issmTime=nstepissm/365; % time where we saved ISSM (yr)
% }}}
disp(['Saving ISSM results to        ' fullfile(resultsrepo,[expname '_issm.mat'])]);
% save ISSM results {{{
issmvars={'issmTime' 'issmVx' 'issmVy' 'issmVel' 'issmThickness' 'issmBasalmeltrate' 'issmX' 'issmY' ...
	'interpVx' 'interpVy' 'interpVel' 'interpThickness' 'interpBasalmeltrate' 'interpX' 'interpY'}; % list of fields 
matname=fullfile(resultsrepo,[expname '_issm.mat']);
save(matname,issmvars{:});
% }}}
% }}}
% gather MITgcm results {{{
% mitgcm coordinates {{{
	LxOC=60E3;  % length of ocean in x (m)
   LyOC=100E3; % length of ocean in y (m)
   Lz=1100;    % depth of ocean in z (m)
   LyICE=60E3; % length of ice domain in y (m)
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
	xg=xp(1:end-1);	% location of x edge points, not including the furthest (redundant) boundary (m)
	yg=yp(1:end-1);   % location of y edge points, not including the furthest (redundant) boundary (m)
   zg=zp(1:end-1);   % location of z edge points, not including the lowest (redundant) boundary (m)

   xc=xp(1:end-1)+0.5*dx; % location of cell centers in x (m)
   yc=yp(1:end-1)+0.5*dy; % location of cell centers in y (m)
   zc=zp(1:end-1)+0.5*dz; % location of cell centers in z (m)

   xcOC=xc(xc>=0 & xc<=LxOC); % location of filled ocean cell centers in x (m)
   ycOC=yc(yc>=0 & yc<=LyOC); % location of filled ocean cell centers in y (m)

   % 2d grid center points
   [XC YC]=meshgrid(xc,yc);
	% 3d grid edge points
   [XG3 YG3 ZG3]=ndgrid(xg,yg,zg);
	% 3d grid center points
	[XC3 YC3 ZC3]=ndgrid(xc,yc,zc);
	indOC=(XC3>=0 & XC3<=LxOC & YC3>=0 & YC3<=LyOC); % XC and YC index for the ocean filled cells
% }}}
% initialize MITgcm data structures {{{
	% 3d fields
	S=zeros(Nx,Ny,Nz,nmitgcm);
	T=zeros(Nx,Ny,Nz,nmitgcm);
	U=zeros(Nx,Ny,Nz,nmitgcm);
	V=zeros(Nx,Ny,Nz,nmitgcm);
	W=zeros(Nx,Ny,Nz,nmitgcm);
	% 2d fields
	Eta=zeros(Nx,Ny,nmitgcm);
	SHICE_fwFluxtave=zeros(Nx,Ny,nmitgcm);
	SHICE_mass=zeros(Nx,Ny,nmitgcm);
% }}}
disp(['Loading MITgcm results from   ' fullfile(expdir,'run')]);
% gather MITgcm results {{{
	y2s=60*60*24*365; % s/yr
	mitgcmdt=str2num(extractAfter(expname,'_dt_')); % mitgcm timestep (s)
	mitgcmTime=nstepmitgcm/365; % time where we saved MITgcm (yr)
	niter=round(mitgcmTime*y2s/mitgcmdt); %mitgcm niter
	rundir=fullfile(expdir,'run'); % where mitgcm results are
	for i=1:nmitgcm
		disp(sprintf('reading step %i/%i',i,nmitgcm));
		niterstr=sprintf('%0.10i',niter(i));
		% 3d fields
		S(:,:,:,i)=binread(fullfile(rundir,['S.' niterstr '.data']),4,Nx,Ny,Nz);
		T(:,:,:,i)=binread(fullfile(rundir,['T.' niterstr '.data']),4,Nx,Ny,Nz);
		% U points (on the Arakawa C grid) are on XG, YC
		% V points (on the Arakawa C grid) are on XC, YG
		Uarakawa=binread(fullfile(rundir,['U.' niterstr '.data']),4,Nx,Ny,Nz);
		F = griddedInterpolant(XG3,YC3,ZC3,Uarakawa);
		U(:,:,:,i)=F(XC3,YC3,ZC3);
		Varakawa=binread(fullfile(rundir,['V.' niterstr '.data']),4,Nx,Ny,Nz);
		F = griddedInterpolant(XC3,YG3,ZC3,Varakawa);
		V(:,:,:,i)=F(XC3,YC3,ZC3);
		W(:,:,:,i)=binread(fullfile(rundir,['W.' niterstr '.data']),4,Nx,Ny,Nz);
		% 2d fields
		Eta(:,:,i)=binread(fullfile(rundir,['Eta.'  niterstr '.data']),4,Nx,Ny);
		SHICE_fwFluxtave(:,:,i)=binread(fullfile(rundir,['SHICE_fwFluxtave.'  niterstr '.data']),4,Nx,Ny);
		SHICE_mass(:,:,i)=binread(fullfile(rundir,['SHICE_mass.'  niterstr '.data']),4,Nx,Ny);
	end
	mitgcmX=XC;
	mitgcmY=YC;
% }}}
disp(['Saving MITgcm results to      ' fullfile(resultsrepo,[expname '_mitgcm.mat'])]);
% save mitgcm results {{{
	mitgcmvars={'mitgcmTime','S','T','U','V','W','Eta','SHICE_fwFluxtave','SHICE_mass','mitgcmX','mitgcmY'};
	matname=fullfile(resultsrepo,[expname '_mitgcm.mat']);
	save(matname,mitgcmvars{:});
% }}}
% }}}
end
