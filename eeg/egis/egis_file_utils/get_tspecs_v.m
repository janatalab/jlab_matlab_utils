function specs = get_tspecs_v(chdr,cell,trial)
%spec_value=get_tspecs_v(cell_hdrs,cell,trial)
%
%Returns a vector which is all the EGIS trial-specifics from cell
%cell_num and trial trial_num. The source  matrix cell_hdrs is
%returned as the second value of the function rd_ses_hdr_v.
%
%See Also:rd_ses_hdr_v, get_tspec_v.
%
%Created 3/95 by Brian Rakitin
%Modified 3/5/96 by BCR
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

numspecs=floor(chdr(cell,LSpec)/2);
tspec_offset = 6+((trial-1)*numspecs); 
specs = chdr(cell,tspec_offset:(tspec_offset+numspecs-1));
