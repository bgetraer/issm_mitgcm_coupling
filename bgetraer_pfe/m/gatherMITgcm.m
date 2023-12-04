function Results=gatherMITgcm(expname,nstep)

% define repositories
exprepo='/nobackupp11/bgetraer/issmjpl/proj-getraer/issm_mitgcm_coupling/bgetraer_pfe/experiments';
expdir=fullfile(exprepo,expname); % where these results are located
rundir=fullfile(expdir,'run'); % where mitgcm results are

disp(['Loading MITgcm results from   ' fullfile(expdir,'run')]);
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
   xg=xp(1:end-1);   % location of x edge points, not including the furthest (redundant) boundary (m)
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
% mitgcm time parameters {{{
y2s=60*60*24*365; % s/yr
mitgcmdt=str2num(extractAfter(expname,'_dt_')); % mitgcm timestep (s)
Results.mitgcmTime=nstep/365; % time where we saved MITgcm (yr)
niter=round(Results.mitgcmTime*y2s/mitgcmdt); %mitgcm niter
% }}}
% gather MITgcm results {{{
disp(sprintf('reading step %i',nstep));
niterstr=sprintf('%0.10i',niter);
% 3d fields
Results.S=binread(fullfile(rundir,['S.' niterstr '.data']),4,Nx,Ny,Nz);
Results.T=binread(fullfile(rundir,['T.' niterstr '.data']),4,Nx,Ny,Nz);
%Results.PH=binread(fullfile(rundir,['PH.' niterstr '.data']),4,Nx,Ny,Nz);
% U points (on the Arakawa C grid) are on XG, YC
% V points (on the Arakawa C grid) are on XC, YG
Uarakawa=binread(fullfile(rundir,['U.' niterstr '.data']),4,Nx,Ny,Nz);
F = griddedInterpolant(XG3,YC3,ZC3,Uarakawa);
Results.U=F(XC3,YC3,ZC3);
Varakawa=binread(fullfile(rundir,['V.' niterstr '.data']),4,Nx,Ny,Nz);
F = griddedInterpolant(XC3,YG3,ZC3,Varakawa);
Results.V=F(XC3,YC3,ZC3);
Results.W=binread(fullfile(rundir,['W.' niterstr '.data']),4,Nx,Ny,Nz);
% 2d fields
Results.Eta=binread(fullfile(rundir,['Eta.'  niterstr '.data']),4,Nx,Ny);
Results.SHICE_fwFluxtave=binread(fullfile(rundir,['SHICE_fwFluxtave.'  niterstr '.data']),4,Nx,Ny);
Results.SHICE_mass=binread(fullfile(rundir,['SHICE_mass.'  niterstr '.data']),4,Nx,Ny);
Results.mitgcmX=xc;
Results.mitgcmY=yc;
Results.mitgcmZ=zc;
% }}}
end
