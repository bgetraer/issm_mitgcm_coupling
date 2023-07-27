function writePARM04(filename,dX,dY,dR,Nx,Ny,Nr,varargin)
%WRITEPARM04 writes CARTESIAN grid options for MITgcm file input/data
% input/data accepts gridding parameters under PARM04. This function assumes an existing input/data file at filename
% and overwrites parameters under PARM04 replacing them with options from the input.
% The following options are set:
%	usingCartesianGrid=.TRUE.
%	delX=Nx*dX
%  delY=Ny*dY
%  delR=Nr*dR
%	IF SPECIFIED:
%	xgOrigin=X0 (default is 0)
%	ygOrigin=y0 (default is 0)
% where Nx*dX is not multiplication, but defines a grid vector with Nx cells with side length of dX.
% Any other parameters defined in PARM04 will throw an error, and this code would need to be updated
% to handle that case.
%
% USE:
%  PARM04(filename,dX,dY,dR,Nx,Ny,Nr);
%  PARM04(filename,dX,dY,dR,Nx,Ny,Nr,'X0',X0,'Y0',Y0);
%  PARM04(filename,dX,dY,dR,Nx,Ny,Nr,'print',1);
% INPUT:
%	filename		location of input/data to be altered
% See also: READPARM04
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
% read until PARM04, write to temp file {{{
   tline = fgetl(readID);
	fprintf(writeID,formatSpec,tline);
	if printfile, disp(tline); end
   while ~contains(tline,'PARM04')
      tline = fgetl(readID);
		fprintf(writeID,formatSpec,tline);
		if printfile, disp(tline); end
   end % }}}
% write new PARM04 to temp file {{{
	% strings to set gridding parameters
	tline={' usingCartesianGrid=.TRUE.,',... % set grid type
		[' delX=' num2str(Nx) '*' num2str(dX) ','],...	% set X grid vector
		[' delY=' num2str(Ny) '*' num2str(dY) ','],...	% set Y grid vector
		[' delR=' num2str(Nr) '*' num2str(dR) ','],...		% set R grid vector
		[' xgOrigin=' num2str(X0) ','],...	% set X origin
		[' ygOrigin=' num2str(Y0) ',']};	% set Y origin
	% write to file
	for i=1:length(tline)
		fprintf(writeID,formatSpec,tline{i});	
		if printfile, disp(tline{i}); end
	end % }}}
	% read until next section of the read file {{{
	tline = fgetl(readID);
	while ~contains(tline,'&')
		tline=tline(~isspace(tline));
		% ignore commented lines
		if tline(1)~='#'
   	   % get the parameter name and value
   	   name=extractBefore(tline,'=');
			if ~any(strcmp({'usingCartesianGrid','delX','delY','delR','xgOrigin','ygOrigin'},name))
				error([filename ' contains PARM04 parameter not supported by WRITEPARM04'])
			end
		end 
	   tline = fgetl(readID);
	end % }}}
% read the rest of the file, write to temp file {{{
while isstr(tline)
	fprintf(writeID,formatSpec,tline);
	if printfile, disp(tline); end
	tline = fgetl(readID);
end % }}}
% }}}
% close both files {{{
   fclose(writeID);
   fclose(readID); % }}}
	% save old file, overwrite filename with the new file {{{
	disp(['Moving existing file ' filename ' to ' filename '_old']);
	eval(['!mv ' filename ' ' filename '_old'])
	disp(['Writing new file ' filename]);
	eval(['!mv ' tmpName ' ' filename]) % }}}
