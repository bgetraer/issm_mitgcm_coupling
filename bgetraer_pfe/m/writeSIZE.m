function writeSIZE(filename,values,varargin)
%WRITESIZE creates a new MITgcm SIZE.h header file from input values.
% Read off of the SIZE.h template provided in MITgcm/model/inc/SIZE.h and write to a new file,
% replacing the variable values with values provided in the input. By default, WRITESIZE saves
% any existing file at filename to filename_old before creating a new file.
% Optional name, value arguments can print the full text of the new SIZE.h file to terminal output
% WRITESIZE to overwrite filename if it already exists.
%	
% USE: 
%	 	writeSIZE(filename,values)
%	 	writeSIZE(filename,values,name,value)
% INPUT:
%	 	values	array of 11 values with the expected order
%	 			sNx, sNy, OLx, OLy, nSx, nSy, nPx, nPy, Nx, Ny, Nr
%	
%	 	varargin
%	 		'print'		bool	print file contents to terminal output
%							use: writeSIZE(values,...,'print',1)
%	 		'overwrite'	bool	overwrite existing file at filepath
%							use: writeSIZE(values,...,'overwrite',1)
% See also READSIZE
% BENJAMIN GETRAER	

	% set default options for print overwrite {{{
	printfile=false; % print file contents to terminal as we go?
	overwritefile=false; % overwrite filepath if it exists? }}}
	% check if values input is useable {{{
	if length(values)~=11
		error('VALUES input must be array of length 11');
	else
		nX=values(1)*values(5)*values(7);
		nY=values(2)*values(6)*values(8);
		if values(9)~=nX || values(10)~=nY
			error('VALUES input inconsistent');
		end
	end % }}}
	% read inputs from varargin if they exist {{{
	for i=1:2:length(varargin)
		switch varargin{i}
			case 'print'
				printfile=varargin{i+1};
			case 'overwrite'
				overwritefile=varargin{i+1};
			otherwise
				error(['Illegal input argument: ' varargin{i}]);
		end
	end % }}}
	% check for existence of file, create backup if necessary {{{
	if exist(filename)==2
		if overwritefile
			disp(['Existing file ' filename ' will be overwritten']);
		else
			disp(['Moving existing file ' filename ' to ' filename '_old']);
			eval(['!mv ' filename ' ' filename '_old']);
		end
	end % }}}
	% open new file and template file {{{
	disp(['Writing new file ' filename])
	writeID=fopen(filename,'w');
	% read from template
	mitgcm_dir='/nobackup/bgetraer/MITgcm';
	template_filename=fullfile(mitgcm_dir, 'model/inc/SIZE.h');
	readID=fopen(template_filename,'r'); % }}}

	% read through the template file, write to the new file {{{
	formatSpec='%s\n'; % new line after each string is written
	% read through any uncommented header, do NOT write to new file {{{
	tline = fgetl(readID);
	while ~strcmp(tline,'CBOP')
		tline = fgetl(readID);
	end % }}}
	% read through the commented header, write to new file {{{
	while tline(1)=='C'
		fprintf(writeID,formatSpec,tline);
		if printfile, disp(tline);end
		tline = fgetl(readID);
	end %}}}
	% read through the variable declarations, write to new file {{{
	while contains(tline,'INTEGER') | contains(tline,'PARAMETER')
		fprintf(writeID,formatSpec,tline);
		if printfile, disp(tline);end
	   tline = fgetl(readID);
	end % }}}
	% read through the variable values, write to new file (assumes default ordering) {{{
	i=1;
	while contains(tline,'&')
		% assumes Nx  = sNx*nSx*nPx and Ny  = sNy*nSy*nPy
		if ~any(i==[9,10]) 
			% extract the string of char before and after the template variable value
			[tempvalue,sline] = regexp(tline,'\d*','Match','split');
			% insert the variable value
			tline=[sline{1} num2str(values(i)) sline{2}];
		end
		% write to new file
		fprintf(writeID,formatSpec,tline);
		if printfile, disp(tline);end
		% advance to the next line
		i=i+1;
	   tline = fgetl(readID);
	end % }}}
% read the rest of the template file, write to new file (assumes MAX_OLX = OLx, and MAX_OLY = OLy) {{{
	while isstr(tline)
		fprintf(writeID,formatSpec,tline);
		if printfile, disp(tline);end
	   tline = fgetl(readID);
	end % }}}
	% }}}
	% close both files {{{
	fclose(writeID);
	fclose(readID); % }}}
