function outchdr = remove_tspecs_v(chdr,cell,trial)
%function outchdr = remove_tspecs_v(chdr,cell,trial)
%
%Returns a matrix which is all the EGIS trial-specifics from chdr
%minus the information from cell and trial. NTrials is decremented
%as well. Intended for use with the matrices returned by rd_egis_hdr_v.
%
%See Also:rd_egis_hdr_v, get_tspecs_v, put_tspecs_v, append_tspecs_v.
%
%Created 5/8/95 by Brian Rakitin
%

%Call the script that defines the EGIS header constants
ses_hdr_offsets_v;

if ~(nargin==3 & nargout==1)
	error('Function get_tspecs_v requires 3 input and 1 output arguments.');
end

if cell > size(chdr,1)
	error(['Cell argument must be no greater than ' int2str(size(chdr,1))]);
end

if trial > chdr(cell,NTrials)
	error(['Trial argument must be no greater than ' int2str(chdr(cell,NTrials))]);
end

outchdr=zeros(size(chdr));
outchdr=chdr;
numspecs=floor(chdr(cell,LSpec)/2);
tspec_offset = 6+((trial-1)*numspecs); 
mask=ones(1,size(chdr,2));
mask(tspec_offset:(tspec_offset+numspecs-1))=zeros(1,size([tspec_offset:(tspec_offset+numspecs-1)],2));
outchdr(cell,1:size(outchdr,2))=zeros(1,size(outchdr,2));
outchdr(cell,1:size(find(mask),2))=chdr(cell,mask);
outchdr(cell,NTrials)=chdr(cell,NTrials)-1;
