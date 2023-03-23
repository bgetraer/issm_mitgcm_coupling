%If needed data files in 
%https://github.com/hgu784/MITgcm_67s/tree/main/initial/input/

%Hard coded parameters
steps=1:3;
clustername='totten';


fbase=[pwd '/'];

interactive=0;
loadonly = 1;
%Cluster parameters{{{
if strcmpi(clustername,'pfe');
	cluster=pfe('grouplist','s1690','cpuspernode',28,'numnodes',1,'time',24*5*60,'interactive',interactive,'processor','bro','queue','long'); %tim in minutes
elseif strcmpi(clustername,'amundsen'),
	cluster=generic('name','amundsen.thayer.dartmouth.edu','np',20,'interactive',0);
else
	disp('cluster not supported yet');
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
if perform(org,'RunSingleCoupleStep'),% {{{

	MITgcmDeltaT=100; % MITgcm time step in seconds
	y2s=31536000; % ye
	rho_ice = 917;

	cd input;
	rdmds_init;
	gendata;
	cd ..

	nPx = 3;
	nPy = 10;
	Nx = 60;
	Ny = 100;
	dx = 1e3;
	nsteps = 5;


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
	time_step = 1/365;
	md = loadmodel(org,'SteadystateNoSlip');

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
        disp(command)
        eval(command)
        newline = [' ntimesteps = ' num2str(time_step*y2s/MITgcmDeltaT)];
        command=['!sed "s/.*ntimesteps.*/' newline '/" data > data.temp; mv data.temp data'];
        disp(command)
        eval(command)
        newline = [' frequency(3) = ' num2str(time_step*y2s)];
        command=['!sed "s/.*frequency(3).*/' newline '/" data.diagnostics > data.temp; mv data.temp data.diagnostics'];
        disp(command)
        eval(command)
        newline = [' pChkptFreq = ' num2str(time_step*y2s)];
        command=['!sed "s/.*pChkptFreq.*/' newline '/" data > data.temp; mv data.temp data'];
        disp(command)
        eval(command)

        results=md.results;

        for coupled_step = 0:(nsteps-1);

         t = time_step * coupled_step;
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


         binwrite([fbase 'run/shelfice_dmdt.bin'],dmdt_mitgcm');
         system(['cp ' fbase 'run/shelfice_dmdt.bin ' fbase 'run/shelfice_dmdt_' num2str(coupled_step) '.bin']);


         newline = [' niter0 = ' num2str(t*y2s/MITgcmDeltaT)];
         command=['!sed "s/.*niter0.*/' newline '/" data > data.temp; mv data.temp data'];
         disp(command)
         eval(command)


         disp('about to run MITgcm')
         tic
         eval(['!mpirun -np ' int2str(nPx*nPy) ' ./mitgcmuv']);
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

         base=md.results.TransientSolution(end).Base;
         thickness=md.results.TransientSolution(end).Thickness;

        end

        savemodel(org,md);
end

