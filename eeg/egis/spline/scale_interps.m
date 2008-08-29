function status=scale_interps(imin,imax,irows,icols,filenums,inroot,outroot,invf);
%function status=scale_interps(imin,imax,irows,icols,filenums,inroot,outroot,invf);
%
%This function takes a 32-bit floating point binary file output
%from interp_by_samp or Interpolator 4.0, rescales it to an 8-bit
%unsigned character, and writes it to a new file.
%

%Created 3/21/96 by B. Rakitin
%

status= -1;

if (nargin~=7) & (nargin~=8)
	error('get_interp_range requires 7 or 8 input params. See help.')
end

n=length(filenums);
if n==0
	error('Null vector passed as filenums param, asshole.');
end

if imin>=imax
	error('Mins must be greater than maxs, dipshit.');
end

if strcmp(inroot,outroot)==1
	error('Input and output paths must be different, you must be Doran.');
end

if nargin==8
	if (invf~=0) & (invf~=1)
		error('Invert flag must be 0 or 1.');
	end
else
	invf=0;
end

tempdata=zeros(irows,icols);

for f=1:n
	infilename=[inroot int2str(filenums(f))];
	disp(['Current File: ' infilename]);
	if exist(infilename)~=2
		disp(['Skipping non-existing input file ' infilename ]);
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
			fclose(infid);
			tempdata=round(((tempdata - imin).*255)/(imax-imin));
			tempdata=tempdata.*(~(tempdata<0));
			tempdata=((tempdata>255)*255)+(tempdata.*(~(tempdata>255)));
%			tempdata(tempdata>255)=ones(size(tempdata(tempdata>255)))*255;
%			tempdata(tempdata<0)=zeros(size(tempdata(tempdata<0)));
			fnumstr=int2str(filenums(f));
			if filenums(f)<100, fnumstr=['0' fnumstr];end;
			if filenums(f)<10, fnumstr=['0' fnumstr];end;
			outfilename=[outroot fnumstr];
			outfid=-1;msg='';
			[outfid,msg]=fopen(outfilename,'w');
			if outfid==-1
				error(['Error opening output file ' outfilename]);
			else
				disp(['Output file is ' outfilename]);
				if invf==0
					count=fwrite(outfid,tempdata,'uchar');
				else
					count=fwrite(outfid,tempdata','uchar');
				end
				if count~=(irows*icols)
					error(['Error writing scaled data to file ' outfilename]);
				end
				fclose(outfid);
			end
		end	%if infid==-1
	end	%if exist(infilename)~=2
end	%for f=1:n

status=1;
