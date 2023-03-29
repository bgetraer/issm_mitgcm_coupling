%If needed data files in 
%https://github.com/hgu784/MITgcm_67s/tree/main/initial/input/

%Hard coded parameters
if ~exist('steps'), steps=1:5; end
clustername='totten';

%parameters
if ~exist('nPx'), nPx=2; end
if ~exist('nPy'), nPy=4; end
Nx = 60;
Ny = 100;
dx = 1e3;
% change as needed
if ~exist('nsteps'); nsteps = 30; end

% change as needed
if ~exist('coupled_time_step'); coupled_time_step = 1/365; end

MITgcmDeltaT=100; % MITgcm time step in seconds
y2s=31536000; % ye
rho_ice = 917;

% the correction parameter: 1 -- fully "corrected", 0 -- no correction

if ~exist('alpha_correction'); alpha_correction = 0; end

fbase=[pwd '/'];

interactive=0;
loadonly = 1;
%Cluster parameters{{{
if strcmpi(clustername,'pfe');
	cluster=pfe('grouplist','s1690','cpuspernode',28,'numnodes',1,'time',24*5*60,'interactive',interactive,'processor','bro','queue','long'); %tim in minutes
elseif strcmpi(clustername,'amundsen'),
	cluster=generic('name','amundsen.thayer.dartmouth.edu','np',20,'interactive',0);
else
%	disp('cluster not supported yet');
end%}}}

%addpath('/local/helene/issm/trunk-jpl/test/MITgcm/tools');
org=organizer('repository',[fbase 'Models'],'prefix','PigLike.','steps',steps);

if perform(org,'GetMITgcm'),% {{{
end%}}}
if perform(org,'Mesh_generation'),% {{{
	md=model();
	md=squaremesh(md,60000,100000,61,101);
	%md=triangle(model(),'./ExpPar/Shelf.exp',1000.);
	md=setmask(md,'all','');
	%Add melt parameterization to run uncoupled
	md=parameterize(md,'./ExpPar/SquareShelf.par');
	md=setflowequation(md,'SSA','all');

	savemodel(org,md);
end%}}}
if perform(org,'SteadystateNoSlip'),% {{{

	md=loadmodel(org,'Mesh_generation');

	md.timestepping.time_step=1/12;
	md.timestepping.start_time=0;
	md.timestepping.final_time=100;
	md.settings.output_frequency=520;

%	md.basalforcings=linearbasalforcings(md);
%	md.basalforcings.deepwater_melting_rate=85.5;
%	md.basalforcings.deepwater_elevation=-750;
%	md.basalforcings.upperwater_melting_rate=0;
%	md.basalforcings.upperwater_elevation=-150;
        pos=find(md.mesh.x<1500 | md.mesh.x==max(md.mesh.x));
	md.stressbalance.spcvy(pos)=0;
	md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.geothermalflux=zeros(md.mesh.numberofvertices,1);
	md.miscellaneous.name='PigLikeNoSlip';

	md.transient.requested_outputs={'default','BasalforcingsFloatingiceMeltingRate','Thickness'};
	md.verbose=verbose('convergence',true,'solution',true,'module',true);
	md=solve(md,'Transient');

	savemodel(org,md);
end%}}}
if perform(org,'RunCouple'),% {{{

	cd input;
	rdmds_init;
	gendata;
	cd ..

	cd run;
	!ln -s ../input/* .
	!ln -s ../build/mitgcmuv .
	!rm data
	!cp ../input/data .
	!rm data.diagnostics
	!cp ../input/data.diagnostics .

	xpoints = 0:dx:Nx*dx;
	ypoints = 0:dx:Ny*dx;
	xpointsmid = .5*(xpoints(1:end-1)+xpoints(2:end));
	ypointsmid = .5*(ypoints(1:end-1)+ypoints(2:end));
	[xpointsmid ypointsmid] = meshgrid(xpointsmid,ypointsmid);
	xpointsmid = xpointsmid(:);
	ypointsmid = ypointsmid(:);
   xpointsmid2 = .5*(xpoints(1:end-1)+xpoints(2:end));
   ypointsmid2 = .5*(ypoints(1:end-1)+ypoints(2:end));

	t=0;
	md = loadmodel(org,'SteadystateNoSlip');
	time_step = coupled_time_step;

	md.results.TransientSolution=md.results.TransientSolution(end);
	base=md.results.TransientSolution(end).Base;
	thickness=md.results.TransientSolution(end).Thickness;
	md.geometry.base=base;
	md.geometry.thickness=thickness;
	md.geometry.surface=md.geometry.base+md.geometry.thickness;
	md.initialization.vx=md.results.TransientSolution(end).Vx;
	md.initialization.vy=md.results.TransientSolution(end).Vy;
	md.initialization.vel=md.results.TransientSolution(end).Vel;
	md.initialization.pressure=md.results.TransientSolution(end).Pressure;

	%md.cluster.executionpath='/local/helene/issmjpl/proj-seroussi/TestCoupling/execution';
	%md.cluster.codepath = '/thayerfs/apps/issm/bin';
	%md.cluster.etcpath = '/thayerfs/apps/issm/etc';
	md.timestepping.final_time=time_step;
	md.timestepping.time_step=1/365;
	md.timestepping.start_time=0;
	md.transient.requested_outputs={'default'};
	md.settings.output_frequency=1;

	md.basalforcings.floatingice_melting_rate(:)=0;

   newline = [' niter0 = ' num2str(t*y2s/MITgcmDeltaT)];
   command=['!sed "s/.*niter0.*/' newline '/" data > data.temp; mv data.temp data'];
   eval(command)
   newline = [' ntimesteps = ' num2str(time_step*y2s/MITgcmDeltaT + 1)];
   command=['!sed "s/.*ntimesteps.*/' newline '/" data > data.temp; mv data.temp data'];
   eval(command)
   newline = [' frequency(3) = ' num2str(time_step*y2s)];
   command=['!sed "s/.*frequency(3).*/' newline '/" data.diagnostics > data.temp; mv data.temp data.diagnostics'];
   eval(command)
   newline = [' frequency(4) = ' num2str(-time_step*y2s)];
   command=['!sed "s/.*frequency(4).*/' newline '/" data.diagnostics > data.temp; mv data.temp data.diagnostics'];
   eval(command)
   newline = [' timephase(4) = ' num2str(0)];
   command=['!sed "s/.*timephase(4).*/' newline '/" data.diagnostics > data.temp; mv data.temp data.diagnostics'];
   eval(command)
   newline = [' pChkptFreq = ' num2str(time_step*y2s)];
   command=['!sed "s/.*pChkptFreq.*/' newline '/" data > data.temp; mv data.temp data'];
   eval(command)

   results=md.results;

   for coupled_step = 0:(nsteps-1);
		md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
		% md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
		% md.basalforcings.geothermalflux=zeros(md.mesh.numberofvertices,1);

		t = time_step * coupled_step;
		if (coupled_step>0);
                 shice_mass_latest = rdmds([fbase 'run/SHICE_mass'], round((t)*y2s/MITgcmDeltaT))';
%                 shice_mass_latest_mesh=InterpFromGridToMesh(xpointsmid2',ypointsmid2',reshape(shice_mass_latest,[Ny,Nx]),md.mesh.x,md.mesh.y,0);
                 issm_mass_latest = md.results.TransientSolution(end).Thickness * rho_ice;
		 issm_mass_latest_grid = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,issm_mass_latest,xpointsmid,ypointsmid,'default',0);
                 dmdt_adjust = alpha_correction * (reshape(issm_mass_latest_grid,[Ny,Nx])-shice_mass_latest) / time_step / y2s;
		 dmdt_adjust(shice_mass_latest==0)=0;
                 
                else
                 dmdt_adjust = zeros(Ny,Nx);
                end

		md.transient.requested_outputs={'default','BasalforcingsFloatingiceMeltingRate','Thickness'};
		tic
		md=solve(md,'Transient');
		toc

		dmdt_icenodes=rho_ice * (md.results.TransientSolution(end).Thickness-md.geometry.thickness)/(md.timestepping.final_time-md.timestepping.start_time);
		% will fail if we include dmdt where there is not meant to be ice
		dmdt_icenodes(md.results.TransientSolution(end).Thickness <=2) = 0;

		% interpolate dmdt from mesh to grid
		dmdt_mitgcm=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,dmdt_icenodes,xpointsmid,ypointsmid,'default',0);
		dmdt_mitgcm=1/y2s*reshape(dmdt_mitgcm,[Ny,Nx]);


		binwrite([fbase 'run/shelfice_dmdt.bin'],dmdt_mitgcm' + dmdt_adjust');
		system(['cp ' fbase 'run/shelfice_dmdt.bin ' fbase 'run/shelfice_dmdt_' num2str(coupled_step) '.bin']);


		newline = [' niter0 = ' num2str(t*y2s/MITgcmDeltaT)];
		command=['!sed "s/.*niter0.*/' newline '/" data > data.temp; mv data.temp data'];
		disp(command)
		eval(command)

		disp('about to run MITgcm')
		tic
		eval(['!mpirun -np ' int2str(nPx*nPy) ' ./mitgcmuv > out 2> err']);
		disp('done MITgcm')
		toc

		system(['cp ' fbase 'run/SHICE_fwFluxtave.' appNum(round((t+time_step)*y2s/MITgcmDeltaT),10) '.data ' fbase 'run/melt.data'])
		system(['cp ' fbase 'run/STDOUT.0000 ' fbase 'run/stdout' num2str(coupled_step)])
		melt=binread('melt.data',4,Nx,Ny)';

		% interpolation of melt onto the issm mesh using xpointsmid and ypointsmid
		melt_mesh=InterpFromGridToMesh(xpointsmid2',ypointsmid2',reshape(melt,[Ny,Nx]),md.mesh.x,md.mesh.y,0);
		md.basalforcings.floatingice_melting_rate=-melt_mesh(:)*y2s/rho_ice;
		md=solve(md,'Transient');

		%Save results of run with melt
		results.TransientSolution(end+1)= md.results.TransientSolution(end);
		results.TransientSolution(end).time = (coupled_step+1)*time_step

		base=md.results.TransientSolution(end).Base;
		thickness=md.results.TransientSolution(end).Thickness;
		md.geometry.base=base;
		md.geometry.thickness=thickness;
		md.geometry.surface=md.geometry.base+md.geometry.thickness;
		md.initialization.vx=md.results.TransientSolution(end).Vx;
		md.initialization.vy=md.results.TransientSolution(end).Vy;
		md.initialization.vel=md.results.TransientSolution(end).Vel;
		md.initialization.pressure=md.results.TransientSolution(end).Pressure;
	end

	md.results = results;
	savemodel(org,md);
	cd ..
end%}}}
if perform(org,'TestDrift'),% {{{ 

	xpoints = 0:dx:Nx*dx;
        ypoints = 0:dx:Ny*dx;
        xpointsmid = .5*(xpoints(1:end-1)+xpoints(2:end));
        ypointsmid = .5*(ypoints(1:end-1)+ypoints(2:end));
        [xpointsmid ypointsmid] = meshgrid(xpointsmid,ypointsmid);
        xpointsmid = xpointsmid(:);
        ypointsmid = ypointsmid(:);
        xpointsmid2 = .5*(xpoints(1:end-1)+xpoints(2:end));
        ypointsmid2 = .5*(ypoints(1:end-1)+ypoints(2:end));


	md = loadmodel(org,'RunCouple');
	for i=0:(nsteps-1);

          if (i>0);
           shelficethick = rdmds('run/SHICE_mass',round(i*y2s*coupled_time_step/MITgcmDeltaT))'/rho_ice;
	  else
           shelficethick = binread('run/shelficemassinit.bin',8,60,100)'/rho_ice;
	  end
	  if (i==0);
	   shelficethick0=shelficethick;
	  end
          issm_thick = md.results.TransientSolution(i+1).Thickness;
          shelficethickmesh = InterpFromGridToMesh(xpointsmid2',ypointsmid2',shelficethick,md.mesh.x,md.mesh.y,0);
	  issm_thick_mitgcm = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,issm_thick,xpointsmid,ypointsmid,'default',0);
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



