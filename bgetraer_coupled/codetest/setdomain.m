% This script controls the MITgcm domain options
% Functions located in issmjpl/proj-getraer/CODE/

% dimensions of the actual filled ocean domain
% filled ocean from 0,LxOC and 0,LyOC and 0,-Lz
LxOC=60E3;	% length of ocean in x (m)
LyOC=100E3;	% length of ocean in y (m)
Lz=1100;		% depth of ocean in z (m)

% resolution in each dimension
dX=1e3;	% m
dY=dX;	% m
dR=20;	% m

% walls in x and y
nWx=2; % number of wall cells in x 
nWy=2; % number of wall cells in y
nWw=ceil(nWx/2);	% number of wall cells before x=0 (excess placed on east)
nWs=nWy;				% number of wall cells before y=0 (excess placed on north)

% origin of full domain
X0=0-nWw*dX;
Y0=0-nWs*dY;

% number of cells in each dimension
NxOC=LxOC/dX;	% number of ocean cells in x
NyOC=LyOC/dY;	% number of ocean cells in y
Nx=NxOC+nWx;	% number of total cells in x
Ny=NyOC+nWy;	% number of total cells in y
Nr=Lz/dR;		% number of total cells in z

% coordinate vectors
xp=((0:Nx)*dX)+X0; % location of all cell edges in x (m)
yp=((0:Ny)*dY)+Y0; % location of all cell edges in y (m)
zp=(0:Nr)*-dR;		 % location of all cell edges in z (m)
xc=xp(1:end-1)+0.5*dX; % location of cell centers in x (m)
yc=yp(1:end-1)+0.5*dY; % location of cell centers in y (m)

% code/SIZE.h parameters
sNx=31;	% Number of X points in tile
sNy=17;	% Number of Y points in tile
OLx=3;	% Tile overlap extent in X
OLy=3;	% Tile overlap extent in Y
nSx=1;	% Number of tiles per process in X
nSy=1;   % Number of tiles per process in Y
nPx=2;	% Number of processes to use in X
nPy=6;   % Number of processes to use in Y

% file locations
fname_SIZE='./code/SIZE.h';
fname_data='./input/data';

% write the SIZE.h file 
values=[sNx, sNy, OLx, OLy, nSx, nSy, nPx, nPy, Nx, Ny, Nr]; % default order (see writeSIZE.m)
writeSIZE(fname_SIZE,values);

% write PARM04 to the input/data file
writePARM04(fname_data,dX,dY,dR,Nx,Ny,Nr,'X0',X0,'Y0',Y0);
