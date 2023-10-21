function draft=scanDepth(thickness)
	rhoi=917; %kg/m^3
	rhow=1028; %kg/m^3
	draft=-thickness*rhoi/rhow;
	mind=min(draft(:));
	maxd=max(draft(:));
	x=0:60;
	for i=1:size(draft,1)
		figure(1); clf;
		plot(x,draft(i,:))
		ylim([mind,maxd]);
		xlim([0,60]);
		title(['y=' num2str(i)])
		pause(0.1);
	end
