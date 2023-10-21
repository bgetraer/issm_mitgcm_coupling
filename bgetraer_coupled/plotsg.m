%md=loadmodel('Models/PigLikeSG.RunCouple');

for i=365*1:20:356*2%length(md.results.TransientSolution)
	H=reshape(md.results.TransientSolution(i).Thickness,61,61);
	V=reshape(md.results.TransientSolution(i).Vel,61,61);
	v=reshape(md.results.TransientSolution(i).Vy,61,61);
	u=reshape(md.results.TransientSolution(i).Vx,61,61);
	M=reshape(md.results.TransientSolution(i).BasalforcingsFloatingiceMeltingRate,61,61);
	D=-H*917/1028;
	h=H(35,:);
	d=-h*917/1028;
	x=0:1E3:60E3;
	%plot(x,d);
	%ylim(-[1000 350]);
	figure(3);clf;
	imagesc(M);colorbar;
	figure(4);clf;
	plot(x,d);
	ylim(-[1000 350]);
	pause(0.05);
end

