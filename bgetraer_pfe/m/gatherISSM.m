function  Results=gatherISSM(expname,nstep)
%GATHERISSM reads results from the file given by expname and nstep and returns a matlab structure
%USEAGE: 
%   Results=gatherISSM(expname,nstep)
% where expname is in the form of 'CTRL_dt_100' for example, and nstep is the coupled step (i.e. 365 for 1 year)
	
exprepo='/nobackupp11/bgetraer/issmjpl/proj-getraer/issm_mitgcm_coupling/bgetraer_pfe/experiments';
expdir=fullfile(exprepo,expname); % where these results are located
disp(['Loading ISSM results from     ' fullfile(expdir,'Models')]);
% interpolation coordinates {{{
dx=1000; %m
Lx=60E3; %m
Ly=60E3; %m
xc=dx/2:dx:Lx;%cell center in x (m)
yc=dx/2:dx:Ly;%cell center in y (m)
[XCice,YCice]=meshgrid(xc,yc);
% }}}
% gather ISSM results {{{
disp(sprintf('reading step %i',nstep));
%load model
filename=sprintf('*%0.5i*Couple.mat',nstep);
filepath=dir(fullfile(expdir,'Models',filename));
md=loadmodel(fullfile(filepath.folder,filepath.name));
%save results on ISSM vertices and elements
Results.issmVx=md.results.TransientSolution(end).Vx; % ice Vx m/yr
Results.issmVy=md.results.TransientSolution(end).Vy; % ice Vy m/yr
Results.issmVel=md.results.TransientSolution(end).Vel; % ice Vel m/yr
Results.issmThickness=md.results.TransientSolution(end).Thickness; % ice thickness m
Results.issmBasalmeltrate=md.results.TransientSolution(end).BasalforcingsFloatingiceMeltingRate; % basal melt rate m/yr
Results.issmX=md.mesh.x; % ice x coordinates (m)
Results.issmY=md.mesh.y; % ice y coordinates (m)
% interpolated onto MITgcm cell center grid
Results.interpVx=griddata(md.mesh.x,md.mesh.y,Results.issmVx,XCice,YCice);
Results.interpVy=griddata(md.mesh.x,md.mesh.y,Results.issmVy,XCice,YCice);
Results.interpVel=griddata(md.mesh.x,md.mesh.y,Results.issmVel,XCice,YCice);
Results.interpThickness=griddata(md.mesh.x,md.mesh.y,Results.issmThickness,XCice,YCice);
elX=mean(md.mesh.x(md.mesh.elements),2); elY=mean(md.mesh.y(md.mesh.elements),2);
Results.interpBasalmeltrate=griddata(elX,elY,Results.issmBasalmeltrate,XCice,YCice);
Results.interpX=XCice;
Results.interpY=YCice;
Results.issmTime=nstep/365; % time where we saved ISSM (yr)
% }}}
end
