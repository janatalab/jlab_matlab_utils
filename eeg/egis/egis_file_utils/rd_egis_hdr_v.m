function [fhdr,chdr,ename,czeros,cgains,cnames,fcom,ftext,coff]=rd_egis_hdr_v(fid)
%[fhdr,chdr,ename,czeros,cgains,cnames,fcom,ftext, coff]=rd_egis_hdr_v(fid)
%
%This reads the header of an EGIS file (either .ave or .ses) into two matrices.
%The first is a row vector containing one element for each field of the 
%description, file info, and session info sections. This vector has room 
%for NCells LCellHdr variables, beggining at column 24. The second matrix
%has NCells rows and a column for each field of the header. The Trial
%Specifics begin in column 6. Room is left for the maximum of LSpec/2 2-byte
%values per trial. The files 'ses_hdr_offsets.m' and 'ave_hdr_offsets.m' provide
%index definitions for the columns of the matrices returned by this function. 
%
%For example the length of the data portion of an EGIS file can be found by
%
%    [my_fhdr, my_chdr]=rd_egis_hdr_v(my)
%    L=my_fhdr(LData)
%
%The ExptNam, Cellname, and flexible area fields are optionally read  
%this function. To retrieve them add return arguments in the following
%order: exp. name, zeros, gains, cell names, comment, text. Arguments
%cannot be skipped, i.e. the first extra argument always returns the
%exp. name.
%
%The function also optionally returns the physical offsets of the beginning of each
%cell's data if coff is passed in the output parameter list
%
%See Also:ses_hdr_offsets_v, ave_hdr_offsets_v, wt_ses_hdr_v, get_tspec_v, get_tspecs_v.

%Modification History:
%	Adapted by Petr Janata from rd_ses_hdr_v.m  6/14/95
%	LHeader read as a 'uint16' 4-30-96 - BCR
%
%Modification History from rd_ses_hdr_v.m
%
%	Orginal Version by Brian Rakitin, 3/2/95
%	Output of rest of header added by BCR, 3/12/95


%Check number of arguments
if nargin ~=1 | fid==-1
	error('Script cut_by_specs requires a valid file id as input.')
end
if nargout == 0 
	error('Script rd_egis_hdr_v requires at least two matrix output arguments.')
else if nargout >= 10
	error('Too many output arguments.')
  end
end

%Prep file for reading
fhdr=zeros(1,23);
status=fseek(fid,0,'bof');
%Read in fhdr fields.
fhdr(1)=fread(fid,1,'int32'); %BytOrd
fhdr(2)=fread(fid,1,'int16'); %HdrVer
fhdr(3)=fread(fid,1,'uint16'); %LHeader
fhdr(4)=fread(fid,1,'int32'); %LData
status=fseek(fid,80,'cof'); %Skip ExptNam
fhdr(5:10)=fread(fid,6,'int16'); %RunDate and RunTime

fhdr(11:23)=fread(fid,13,'int16'); %Everything else up to LCellHdr

fhdr(24:(24+fhdr(18)-1))=fread(fid,fhdr(18),'int16'); %LCellHdr for each cell

%Read in chdr
if nargout >=2 
	chdr=zeros(fhdr(18),5);
	status=fseek(fid,(4*fhdr(19)),'cof');
	for loop=1:fhdr(18)
		chdr(loop,1)=fread(fid,1,'int16');
		status=fseek(fid,80,'cof');
		chdr(loop, 2:5)=(fread(fid,4,'int16'))';
		lspectot=(chdr(loop,2)*(chdr(loop,5)/2));
		lastcol=(6+lspectot-1);
		if lastcol >= 6
		  chdr(loop,lastcol)=0;
		  chdr(loop,6:lastcol)=(fread(fid,lspectot,'int16'))';
		end
	end
else
	return
end

%Read experiment name if necessary
if nargout >= 3 
	status=fseek(fid,12,'bof');
	tempstr=fread(fid,80,'char');
	ename=setstr(tempstr);
else 
	return
end

%Read  zeros if necessary
if nargout >= 4
	status=fseek(fid,(130+(2*fhdr(18))),'bof');
	czeros=fread(fid,(fhdr(19)),'int16');
else
	return
end

%Read gains if necessary
if nargout >= 5
	status=fseek(fid,(130+(2*fhdr(18))+(2*fhdr(19))),'bof');
	cgains=fread(fid,(fhdr(19)), 'int16');
else
	return
end

%Read cellnames if necesaary
if nargout >= 6
	status=fseek(fid,(130+(2*fhdr(18))+(4*fhdr(19))),'bof');
	status=fseek(fid,2,'cof');
	for loop=1:fhdr(18)
		cnames(loop,:)=setstr(fread(fid,80,'char')');
		status=fseek(fid,fhdr(24+(loop-1))-80,'cof');
	end
else
	return
end
		
%Read comment if necessary
if nargout >= 7
	status=fseek(fid,fhdr(3),'bof');
	status=fseek(fid,-(fhdr(22)+fhdr(21)+fhdr(20)),'cof');
	fcom=fread(fid,fhdr(20),'char');
else
	return
end

%Read text if necessary
if nargout >= 8
	status=fseek(fid,fhdr(3),'bof');
	status=fseek(fid,-(fhdr(22)+fhdr(21)),'cof');
	ftext=fread(fid,fhdr(21),'char');
else
	return
end

%Read cell offsets if necessary
if nargout >= 9
	coff=get_cell_offsets(fhdr, chdr);
else
	return
end

