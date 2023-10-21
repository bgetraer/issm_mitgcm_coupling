function draft=thick2draft(thickness)
	rho_i=917;	%kg/m^3
	rho_w=1028;	%kg/m^3
	draft=-thickness*rho_i/rho_w;
