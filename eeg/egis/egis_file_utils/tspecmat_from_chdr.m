function tsmat = tspecmat_from_chdr(chdr)
% tsmat = tspecmat_from_chdr(chdr)
%
% Creates a matrix of tspec values
%
% Format of output matrix:
%
% Column 1: Cell ID
% Column 2: Entry number
% Column 3: Edit code
% 

% 01/09/03  PJ  Fixed reshape command to be able to handle heterogeneous
%               numbers of trials.

ses_hdr_offsets_v;

ncells = size(chdr,1);

ntrials = chdr(:,NTrials);
nspec = chdr(:,LSpec)/2;

max_nspec = max(nspec);

tsmat = [];

for icell = 1:ncells
  tmp = reshape(chdr(icell,TSpecs:(TSpecs+max_nspec*ntrials(icell)-1)),max_nspec,ntrials(icell))';
  
  tsmat = [tsmat; ones(ntrials(icell),1)*icell (1:ntrials(icell))' tmp];
end % for icell=