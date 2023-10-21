%md=loadmodel('../coupled_experiment/Models/PigLike.SteadystateNoSlip');
md=loadmodel('./Models/PigLike.SteadystateNoSlip');
y=md.mesh.y;
posy=(md.mesh.x==30E3);
x=md.mesh.x;
posx=(md.mesh.y==60E3);

figure(1);clf;hold on;
figure(2);clf;hold on;
for i=1:length(md.results.TransientSolution)
	H=md.results.TransientSolution(i).Thickness;
	figure(1);
	plot(y(posy),H(posy));

	figure(2);
	plot(x(posx),H(posx));
end
