function PARM04=readPARM04(filename,varargin)
%READPARM04 reads CARTESIAN grid options for MITgcm file input/data.
% input/data accepts gridding parameters under PARM04. This function reads the following parameters
% declared under PARM04 and returns them as a structure:
%		usingCartesianGrid
%		delX
%		delY
%		delR
%		xgOrigin
%		ygOrigin
% Any other parameters defined in PARM04 will throw an error, and this code would need to be updated
% to handle that case.
%
% USE:
%	PARM04=writePARM04(filename);
%  PARM04=writePARM04(filename,name,value);
% INPUT:
%	filename			location of input/data to be read
%	varargin
%		'print'	0,1	option to print the output of PARM04 to terminal (for debugging). default is 0.
% OUTPUT:
%	PARM04		structure containing the following fields:
%		usingCartesianGrid	(this should always equal ".TRUE.")
%		delX	str	declaration for the X grid vector	
%		delY	str	declaration for the Y grid vector	
%		delR	str	declaration for the Z grid vector
%		dX		num	cell length in X
%		dY		num	cell length in Y
%	`	dR		num	cell length in Z
%		Nx		num	number of cells in X
%     Ny    num   number of cells in Y
%     Nr    num   number of cells in Z
% See also: WRITEPARM04
%
% Benjamin Getraer

% check if print contents to terminal {{{
   printfile=0;
   for i=1:2:length(varargin)
      switch varargin{i}
         case 'print'
            printfile=varargin{i+1};
         otherwise
            error(['Illegal input argument: ' varargin{i}]);
      end
   end % }}}
% open file for reading {{{
	% open file
	fileID=fopen(filename,'r'); % }}}
% read through until PARM04 {{{
   tline = fgetl(fileID);
   while ~contains(tline,'PARM04')
      tline = fgetl(fileID);
		%disp(tline)
   end % }}}
% read PARM04, save values to output structure {{{
if printfile, disp('Reading PARM04:'); end
tline = fgetl(fileID);
while ~contains(tline,'&')
	% eliminate white space
	tline=tline(~isspace(tline));
	% ignore commented lines
	if tline(1)~='#'
		% get the parameter name and value
		name=extractBefore(tline,'=');
		switch name
			case 'usingCartesianGrid'
				PARM04.usingCartesianGrid=string(extractBetween(tline,'=',','));
			case 'delX'
				PARM04.delX=string(extractBetween(tline,'=',','));
				PARM04.Nx=str2num(string(extractBefore(PARM04.delX,'*')));
				PARM04.dX=str2num(string(extractAfter(PARM04.delX,'*')));
			case 'delY'
				PARM04.delY=string(extractBetween(tline,'=',','));
				PARM04.Ny=str2num(string(extractBefore(PARM04.delY,'*')));
            PARM04.dY=str2num(string(extractAfter(PARM04.delY,'*')));
			case 'delR'
				PARM04.delR=string(extractBetween(tline,'=',','));
				if contains(PARM04.delR,'*')
					PARM04.Nr=str2num(string(extractBefore(PARM04.delR,'*')));
					PARM04.dR=str2num(string(extractAfter(PARM04.delR,'*')));
				else
					PARM04.Nr=1;
					PARM04.dR=str2num(PARM04.delR);
				end
			case 'xgOrigin'
				PARM04.xgOrigin=string(extractBetween(tline,'=',','));
			case 'ygOrigin'
				PARM04.ygOrigin=string(extractBetween(tline,'=',','));
			otherwise
				error([filename ' contains PARM04 parameter not supported by READPARM04'])
		end
		if printfile, disp(tline); end
	end
	tline = fgetl(fileID);
end % }}}
