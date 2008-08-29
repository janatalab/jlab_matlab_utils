function outchdr= append_tspecs_v(inspecs,chdr,cell)
%outchdr=append_tspecs_v(inspecs,chdr,cell)
%
%For use with the EGIS file utilities. 
%
%Writes the vector inspecs to the structure chdr at the 
%position specified by cell. Assumes the trial is NTrials
%+1, and updates this value. Expands the chdr structure if
%necessary.  
%
%See Also:rd_ses_hdr_v, get_tspec_v.
%
%Created 3/5/96 by Brian Rakitin
%Udating of NTrials in chdr fixed 5/4/96 BCR
%

%Call the script that defines the EGIS header constants
ses_hdr_offsets_v;

if ~(nargin==3 & nargout==1)
	error('Function append_tspecs_v requires 3 input and 1 output arguments.');
end

if cell==0
	error('Cell argument cannot be zero.');
end

if cell > size(chdr,1)
	error(['Cell argument must be no greater than ' int2str(size(chdr,1))]);
end

if size(inspecs)~=floor(chdr(cell,LSpec)/2)
	error(['Input argument inspecs must be length ' int2str(floor(chdr(cell,LSpec)/2))]);
end

trial=chdr(cell,NTrials)+1;
chdr(cell,NTrials)=trial;
numspecs=floor(chdr(cell,LSpec)/2);
tspec_offset = 6+((trial-1)*numspecs); 
if (tspec_offset+numspecs-1) <= size(chdr,2)
	outchdr=chdr;
else
	outchdr=zeros(size(chdr,1),(trial*numspecs)+5);
	outchdr(1:size(chdr,1),1:size(chdr,2))=chdr;
end
outchdr(cell,tspec_offset:(tspec_offset+numspecs-1))=inspecs;
