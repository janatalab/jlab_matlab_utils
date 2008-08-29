function [intmin, intmax]=get_interp_range(irows,icols,filenums,inroot);
%function [intmin, intmax]=get_interp_range(irows,icols,filenums,inroot);
%
%Parameters:
%	irows	- Total rows per image file
%	icols	- Total cols per image file
%	filenums	- vector of numbers appended to inroot to create
%		input file names
%	inroot	- Used with above to create input file names
%
%This function scales images created with interp_by_samp.m using
%the 's' header_flag parameter. Goes through the files specified by
%a combination of the inroot with a filenumber appended, and returns
%the minimum and maximum values found in those files in the output 
%parameters.
%

%Created 3/20/96 by B. Rakitin
%

%disp(['Nargin = ' int2str(nargin) ', Nargout = ' int2str(nargout)]);

if ~(nargin==4)
	error('get_interp_range requires four input params. See help.')
end

if ~((nargout==1)|(nargout==2))
	error('get_interp_range requires one or two output params. See help.')
end

n=length(filenums);
if n==0
	error('Null vector passed as filenums param, asshole.');
end

tempdata=zeros((irows * icols),1);

for f=1:n
	infilename=[inroot int2str(filenums(f))];
	disp(['Current File: ' infilename]);
	if exist(infilename)~=2
		disp(['Skipping non-existing file ' infilename ]);
	else
		infid=-1;msg='';
		[infid,msg]=fopen(infilename);
		if infid==-1
			disp([msg ': ' infilename]);
		else
			[tempdata, count]=fread(infid,size(tempdata),'float');
			if count~=(irows*icols)
				fclose('all');
				error(['Error reading file ' infilename]);
			end
			if f==1
				disp('Initializing temp values.')
				tempmax=max(tempdata);
				tempmin=min(tempdata);
			else
				if max(tempdata)>tempmax,tempmax=max(tempdata);end
				if min(tempdata)<tempmin,tempmin=min(tempdata);end
			end
			fclose(infid);
			disp(['Current Max = ' int2str(tempmax) ', Current Min = ' int2str(tempmin)]);
		end	%if infid==-1
	end	%if exist(infilename)~=2
end	%for f=1:n

if nargout==2
	intmax=tempmax;
	intmin=tempmin;
else
	if nargout==1
		intmin=zeros(1,2);
		intmin(1)=tempmax;
		intmin(2)=tempmin;
	end
end


