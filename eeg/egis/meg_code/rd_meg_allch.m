function trialdata = rd_meg_allch(fid,header_length,NChan,numsamp,start_samp);
%trialdata = rd_meg_allch(fid,header_length,NChan,numsamp,start_samp);
%meg
%reads in meg data from .meg file.  Can be used in two ways 
% one way is to read from an arbitrary start point in samples, start_samp
%	a total of numsamp samples. 
%       in this case file_position is set to 'bof' on calls 
% the other way is to sequentially read the blocks in which case the 
%  read starts at the beginning of the file in blocks of numsamp 
%	in thsi case file_position i to 'cof';
%fid must be a valid header; 
% 
%

if ~(nargin == 4|nargin == 5);
	error('improper number of input arguments')
end;


if nargin == 5
	frewind(fid);
	fseek(fid,header_length,'bof');
	skip_bytes = NChan*(start_samp-1)*4;
	fseek(fid,skip_bytes,'cof');
end;
trialdata  = fread(fid,[NChan,numsamp],'real*4');
trialdata = trialdata';                                       





