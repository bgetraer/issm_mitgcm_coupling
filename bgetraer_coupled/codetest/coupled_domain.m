domain=struct;
% MITGCM DOMAIN
domain.mit.dx=1E3; % horizontal mesh resolution (km)
domain.mit.nx_ocean=60E3/domain.mit.dx; % number of ocean domain cells in x
domain.mit.ny_ocean=100E3/domain.mit.dx; % number of ocean domain cells in y
domain.mit.nx=domain.mit.nx_ocean+2; % total number of cells in x
domain.mit.ny=domain.mit.ny_ocean+2; % total number of cells in y
% throw error if total domain smaller than ocean domain
if (domain.mit.nx<domain.mit.nx_ocean || domain.mit.ny<domain.mit.ny_ocean)
    error('Domain must be at least as large as the ocean domain.')
end

%divide x wall cells between west and east sides, placing first on west
domain.mit.x0=ceil((domain.mit.nx-domain.mit.nx_ocean)/2); % number of wall cells before x=0
%place y wall cells at the southern side 
domain.mit.y0=domain.mit.ny-domain.mit.ny_ocean; % number of wall cells before y=0
domain.mit.x=((0:domain.mit.nx)-domain.mit.x0)*domain.mit.dx; % location of cells in x (km)
domain.mit.y=((0:domain.mit.ny)-domain.mit.y0)*domain.mit.dx; % location of cells in y (km)

domain.mit.outerdomain=[([0,0,domain.mit.nx,domain.mit.nx]-domain.mit.x0)' ...
    ([0,domain.mit.ny,domain.mit.ny,0]-domain.mit.y0)']*domain.mit.dx; % outline of full domain
domain.mit.oceandomain=[[0,0,domain.mit.nx_ocean,domain.mit.nx_ocean]' ...
    [0,domain.mit.ny_ocean,domain.mit.ny_ocean,0]']*domain.mit.dx; % outline of ocean domain
domain.mit.walldomain = [domain.mit.outerdomain;nan,nan;domain.mit.oceandomain]; % outline of walls
% ISSM DOMAIN
domain.issm.dx=domain.mit.dx; % horizontal mesh resolution (km)
domain.issm.nx=domain.mit.nx_ocean; % total number of square elements in x
domain.issm.ny=60E3/domain.issm.dx; % total number of square elements in y
% throw error if ice domain larger than ocean domain
if (domain.issm.ny>domain.mit.ny_ocean)
    error('Ice must be smaller than or equal to the ocean domain.')
end

domain.issm.x=(0:domain.issm.nx)*domain.issm.dx; % location of cells in x (km)
domain.issm.y=(0:domain.issm.ny)*domain.issm.dx; % location of cells in y (km)

domain.issm.icedomain=[([0,0,domain.issm.nx,domain.issm.nx])' ...
    ([0,domain.issm.ny,domain.issm.ny,0])']*domain.issm.dx; % outline of ice domain


% plot the domain
figure(1);clf;hold on;
h_wall = plot(polyshape(domain.mit.walldomain));
h_ice = plot(polyshape(domain.issm.icedomain));
xline(domain.mit.x);
yline(domain.mit.y);
axis equal
axis([min(domain.mit.x) max(domain.mit.x) min(domain.mit.y),max(domain.mit.y)])
set(gca,'XTickLabel',get(gca,'XTick').*1E-3,...
    'YTickLabel',get(gca,'YTick').*1E-3)
xlabel('X (km)');ylabel('Y (km)');
legend([h_wall,h_ice],'Wall','Ice')
title(sprintf('%i by %i cell domain',domain.mit.ny,domain.mit.nx))
%exportgraphics(gcf,'~/Desktop/domain.pdf')
