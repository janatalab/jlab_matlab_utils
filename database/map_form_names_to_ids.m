function [ids] = map_form_names_to_ids(names)
% Calls form_name_defs an maps form names to their IDs
%
% [ids] = map_form_names_to_ids(names,map);
%

form_name_defs;  % Get list of form name mappings

if ~iscell(names)
  names = {names};
end

nid = length(names);
ids = zeros(nid,1);

% Get a listing of the field names 
[form_list_mask, form_list_loc] = ...
    ismember(names,{form_name_id_const_map{:,1}});

for iid = 1:nid
  loc = form_list_loc(iid);
  if loc
    ids(iid) = form_id_const.(form_name_id_const_map{loc,2});
  end
end % for iid
