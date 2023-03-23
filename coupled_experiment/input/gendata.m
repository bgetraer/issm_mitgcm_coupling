%addpath /home/toshiki/MITgcm/utils/matlab/cs_grid/read_cs/
%addpath /home/toshiki/MITgcm/utils/matlab/cs_grid/
%addpath /home/toshiki/MITgcm/utils/matlab/
%%
% cd ~/coupled/input/

% need to chang this part according to the ice-only run you want to use as
% an initial thickness

% change and uncomment
% surf = rdmds('~/initial/run/surfDiag.0000001000'); % <-- change 
% h = surf(:,:,3);
% fid = fopen('H_new.box','w','b'); fwrite(fid,h,'real*8'); fclose(fid);

% for now, i have tested this with the H_new.box in the folder

function gendata

rho_ice = 917;

%% Dimensions of grid
nx=60; 
ny=100;
nz=55;
delz = 20;
dx = 1e3;


hfacMin = 0.2;
%mwct = 3;


%eos = 'linear';
eos = 'jmd95z';
% eos = 'mdjwf';

acc = 'real*8';

dz = delz*ones(1,nz);
zgp1 = [0,cumsum(dz)];
zc = .5*(zgp1(1:end-1)+zgp1(2:end));
zg = zgp1(1:end-1);
dz = diff(zgp1);

% xgOrigin = -105.5,
% ygOrigin =  -75.4457,
% delX=60*.03125,
% delY=100*.0078125,

lat0 = -75.5557;
lon0 = -105.5;
lat = lat0:.0078125:(lat0+100*.0078125);
latc = .5*(lat(1:end-1)+lat(2:end));
lon = lon0:.03125:(lon0+60*.03125);
lonc = .5*(lon(1:end-1)+lon(2:end));
[lon lat]=meshgrid(lon,lat);
lon = lon(:);
lat = lat(:);
save latlon.mat lon lat

xpoints = 0:dx:nx*dx;
ypoints = 0:dx:ny*dx;
xpointsmid = .5*(xpoints(1:end-1)+xpoints(2:end));
ypointsmid = .5*(ypoints(1:end-1)+ypoints(2:end));
[xpointsmid ypointsmid] = meshgrid(xpointsmid,ypointsmid);
xpointsmid = xpointsmid(:);
ypointsmid = ypointsmid(:);


%%% HELENE: HERE WE NEED TO REPLACE WITH READING A BINARY FILE FROM THE INITIAL THICKNESS OF STREAMICE
%%% AND THIS NEEDS TO GO INTO VARIABLE VV %%%%%%%%%%%%%%%%%%%%%%a

md = loadmodel('../Models/PigLike.SteadystateNoSlip');
thickness_mitgcm=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(end).Thickness,xpointsmid,ypointsmid,'default',0);
VV = reshape(thickness_mitgcm,[ny nx])';
binwrite('H_new.box',VV);
%fid = fopen('../codeinput_piglike_streamice/input/H_new.box','r','b'); VV2 = fread(fid,[nx ny],'real*8'); fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bathy = -1100 * ones(ny,nx); 
bathy(1,:) = 0; bathy(:,[1 end]) = 0;
fid = fopen('BATHY60.box','w','b'); fwrite(fid,bathy','real*8'); fclose(fid);

%%%%% SO WE CAN SEE BOUNDARY CONDITIONS FROM STREAMICE
%avgV = 1000;
%x = dx/2:dx:(nx*dx-dx/2);
%xmid = nx*dx/2;
%Vin = (2*avgV/xmid^2) * (xmid^2-(x - xmid).^2);
%vdirich = zeros(ny,nx);
%vdirich(2,:) = Vin;
%fid = fopen('VDIRICH60.box','w','b'); fwrite(fid,vdirich','real*8'); fclose(fid);

%x = 0 : (5*pi)/59 : 5*pi; hin = 1000 + 50*sin(x);
%y = 0 : (13*pi)/59 : 13*pi ; Hin = hin + (150*sin(y));
% Hin = 1200;
%Hbc = zeros(ny,nx);
%Hbc(2,:) = Hin;
%fid = fopen('HBCy60.box','w','b'); fwrite(fid,Hbc','real*8'); fclose(fid);


%%%%%%%%% stratification %%%%%%%%%%%%%%%%

fid = fopen('theta.init','r','b'); q = fread(fid,inf,'real*8'); fclose(fid); Tinit=reshape(q,[nx ny nz]);
fid = fopen('salt.init','r','b'); q = fread(fid,inf,'real*8'); fclose(fid); Sinit=reshape(q,[nx ny nz]);
sref = squeeze(Sinit(1,1,:));
tref = squeeze(Tinit(1,1,:));

%%%%%%%%%%% density %%%%%%%%%%%%%%%%%

% Gravity
gravity=9.81;
rhoConst = 1000;
k=1;
dzm = abs([zg(1)-zc(1) .5*diff(zc)]);
dzp = abs([.5*diff(zc) zc(end)-zg(end)]);
p = abs(zc)*gravity*rhoConst*1e-4;
dp = p;
kp = 0;

Rho = zeros(nz,1);

while rms(dp) > 1e-13
  phiHydF(k) = 0;
  p0 = p;
  kp = kp+1;
  for k = 1:nz
    switch eos
     case 'linear'

     case 'jmd95z'
      drho = densjmd95(sref(k),tref(k),p(k))-rhoConst;
     case 'mdjwf'
      drho = densmdjwf(sref(k),tref(k),p(k))-rhoConst;
     otherwise
      error(sprintf('unknown EOS: %s',eos))
    end
    Rho(k) = drho+rhoConst;
    phiHydC(k)   = phiHydF(k) + dzm(k)*gravity*drho/rhoConst;
    phiHydF(k+1) = phiHydC(k) + dzp(k)*gravity*drho/rhoConst;
  end
  switch eos
   case 'mdjwf'
    p = (gravity*rhoConst*abs(zc) + phiHydC*rhoConst)/gravity/rhoConst;
  end
  dp = p-p0;
end

%shelficemass = binread('HINIT60.bin',nx,ny) * 917;


VV(1,:) = 0;
VV(:,1) = 0;
shelficemass=VV*rho_ice;

topo = zeros(nx,ny);


for ix=1:nx
  for iy=1:ny

     mass = shelficemass (ix,iy);
     massFuncC = rhoConst * (phiHydC/gravity + zc);
     massFuncF = rhoConst * (phiHydF/gravity + zgp1);

     k = max (find ( massFuncF < mass ));
     
     if (isempty(k))
         k=0;
     end
     if (k>0)
     if (mass < massFuncC(k))
      topo(ix,iy) = -zg(k) - (mass-massFuncF(k)) * delz/2 / (massFuncC(k)-massFuncF(k));
     else
      topo(ix,iy) = -zc(k) - (mass-massFuncC(k)) * delz/2 / (massFuncF(k+1)-massFuncC(k));
     end
     end

  end
end

etainit = zeros(size(topo));

% new topography: icetopo rounded to the nearest k * deltaZ
%                 eta_init set to make difference

icetopo2 = topo;

for ix=1:nx
  for iy=1:ny
    k=max(find(abs(zg)<abs(icetopo2(ix,iy))));
    if isempty(k)
      k=0;
    else
      
      dr = 1-(-zg(k) - icetopo2(ix,iy))/delz;
      if (dr > .25)
          % bring Ro_surf *up* to closest grid face & make etainit negative
          % to compensate
          icetopo2(ix,iy) = -zg(k);
          etainit(ix,iy) = (dr-1)*delz;
      else
          % bring Ro_surf *down* to closest grid face & make etainit pos
          % to compensate
          icetopo2(ix,iy) = -zg(k+1);
          etainit(ix,iy) = (dr)*delz;
      end
       
    end
  end
end

%%% runoff
flux = zeros(size(shelficemass));
flux(30,2) = 200.;
 

%etainit(:,:)=0;
%icetopo2(:,1)=0;
fid = fopen('shelftopo.round.bin','w','b'); fwrite(fid,icetopo2,'real*8'); fclose(fid);
fid = fopen('etainit.round.bin','w','b'); fwrite(fid,etainit,'real*8'); fclose(fid);
fid = fopen('shelficemassinit.bin','w','b'); fwrite(fid,shelficemass,'real*8'); fclose(fid);
fid = fopen('runoff_flux.bin','w','b'); fwrite(fid,flux,'real*8'); fclose(fid);
%fid = fopen('bathy_step.bin','w','b'); fwrite(fid,bathy,'real*8'); fclose(fid);
return
