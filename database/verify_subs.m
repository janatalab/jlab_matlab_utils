function expinfo = verify_subs(expinfo,proc_subs)

% Removes subjects associated with the experiment info who aren't on the list
% in proc_subs
%
% expinfo = verify_subs(expinfo,proc_subs);
%

% 10/29/06 Petr Janata

nexp = length(expinfo);

for iexp = 1:nexp
  bad_sub_idxs = find(~ismember(expinfo(iexp).subs.ids,proc_subs));
  
  if ~isempty(bad_sub_idxs)
    nbad = length(bad_sub_idxs);
    
    fprintf('\tRemoving data for %d subjects not on inclusion list\n', nbad);
    
    expinfo(iexp).subs.names(bad_sub_idxs) = [];
    expinfo(iexp).subs.ids(bad_sub_idxs) = [];
    expinfo(iexp).subs.nsess(bad_sub_idxs) = [];
    expinfo(iexp).subs.ages(bad_sub_idxs) = [];
    
    if isfield(expinfo.subs,'sess')
      expinfo(iexp).subs.sess(bad_sub_idxs) = [];
    end
    
  end % if ~isempty(bad_sub_idxs)
end % for iexp=
