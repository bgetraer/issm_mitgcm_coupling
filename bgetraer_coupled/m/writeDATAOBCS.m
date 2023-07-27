function writeDATAOBCS(filename,Nx,Ny,varargin)
%WRITEDATAOBCS edits the domain dependent variables for MITgcm file input/data.obcs
% input/data.obcs accepts parameters for the open boundary condition package. This function
% reads in an existing data.obcs file and overwrites the OB_Jnorth parameter to adjust the 
% boundary based on the size of the domain.
% The parameter is rewritten as follows:
%	OB_Jnorth=Nx*Ny
%
% USE:
%	writeDATAOBCS('data.obcs',Nx,Ny)
% INPUT:
%	filename		location of input/data.obcs to be altered
%	Nx				number of cells in x dimension
%	Ny				number of cells in y dimension
%	varargin
%		'print',1	print out the new file contents as it is written
%
% Benjamin Getraer%

% read varargin {{{
	% set defaults
	X0=0;
	Y0=0;
   printfile=0;
   for i=1:2:length(varargin)
      switch varargin{i}
			case 'X0'
				X0=varargin{i+1};
			case 'Y0'
            Y0=varargin{i+1};
         case 'print'
            printfile=varargin{i+1};
         otherwise
            error(['Illegal input argument: ' varargin{i}]);
      end
   end % }}}
% open files for reading and writing {{{
   % open the existing input/data file
   readID=fopen(filename,'r'); 
	% open a temporary file to write to
	tmpName = tempname;
	writeID=fopen(tmpName,'w'); % }}}
% read through existing file, write to temp file {{{
   formatSpec='%s\n'; % new line after each string is written
   tline = fgetl(readID);
   while ~contains(tline,'OB_Jnorth')
		fprintf(writeID,formatSpec,tline);
		if printfile, disp(tline); end
      tline = fgetl(readID);
   end 
	% write new OB_Jnorth
	tline=[' OB_Jnorth=' num2str(Nx) '*' num2str(Ny) ','];
	fprintf(writeID,formatSpec,tline);
	if printfile, disp(tline); end
	% read the rest of the file, write to temp file
	tline = fgetl(readID);
	while isstr(tline)
		fprintf(writeID,formatSpec,tline);
		if printfile, disp(tline); end
		tline = fgetl(readID);
	end
% }}}
% close both files {{{
   fclose(writeID);
   fclose(readID); % }}}
	% overwrite filename with the new file {{{
	disp(['Writing file ' filename]);
	eval(['!mv ' tmpName ' ' filename]) % }}}
