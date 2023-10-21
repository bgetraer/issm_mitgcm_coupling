function SZ=readSIZE(filename,varargin)
%READSIZE reads an MITgcm SIZE.h header file and outputs the variables to a struct SZ.
%  Optional name, value argument prints the text of SIZE.h to the terminal output as it is read.
%
%  USE:
%     SZ=readSIZE(filename)
%     SZ=readSIZE(filename,name,value)
%  INPUT:
%     filename			location of SIZE.h file to be read      
%     varargin
%        'print'        bool  print file contents to terminal output
%                             use: writeSIZE(values,...,'print',1)	
%
% OUTPUT:
%	SZ			structure with two fields, showing the name and values of the SIZE.h variables
%		SZ.name			sNx, sNy, OLx, OLy, nSx, nSy, nPx, nPy, Nx, Ny, Nr
%		SZ.values		an array of length 11 of values corresponding to SZ.name
%
%  See also WRITESIZE
%  BENJAMIN GETRAER
	
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
	% initialize output structure and open file for reading {{{
	SZ=struct;	
	SZ.name={'sNx','sNy','OLx','OLy','nSx','nSy','nPx','nPy','Nx','Ny','Nr'};
	SZ.values=[];
	% open file
	fileID=fopen(filename,'r'); % }}}
	% read through any uncommented header {{{
	tline = fgetl(fileID);
	while ~strcmp(tline,'CBOP')
		tline = fgetl(fileID)
	end % }}}
	% read through the commented header {{{
	while tline(1)=='C'
		if printfile, disp(tline);end
		tline = fgetl(fileID);
	end % }}}
	% read through the variable declarations {{{
	while contains(tline,'INTEGER') | contains(tline,'PARAMETER')
		if printfile, disp(tline);end
	   tline = fgetl(fileID);
	end % }}}
	
	% read the variable values (assumes default ordering) {{{
	i=1;
	while contains(tline,'&')
		% extract the value of the variable as well as the string of char before and after
		[value,sline] = regexp(tline,'\d*','Match','split');
		% match the expected variable name to the string we read
		if ~contains(sline{1},SZ.name{i})
			error('variables read from SIZE.h are not in expected order.')
		end
		% save the values
		if i==9
			% assumes Nx  = sNx*nSx*nPx
			SZ.values(i)=SZ.values(1)*SZ.values(5)*SZ.values(7);
		elseif i==10
			% assumes Ny  = sNy*nSy*nPy
	      SZ.values(i)=SZ.values(2)*SZ.values(6)*SZ.values(8);
		else
			% save the value
			SZ.values(i)=str2num(value{1});
		end
		% advance to the next line
		i=i+1;
		if printfile, disp(tline);end
	   tline = fgetl(fileID);
	end % }}}
	% read the rest of the file (assumes MAX_OLX = OLx, and MAX_OLY = OLy, does not save them) {{{
	while isstr(tline)
		if printfile, disp(tline);end
	   tline = fgetl(fileID);
	end % }}}
	% close file {{{
	fclose(fileID); % }}}
