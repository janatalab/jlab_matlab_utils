function s = id_map2struct(id_map)
% Creates a struct field/value pairs from the first/second column of the input. 
% s = id_map2struct(id_map);
%
% Generates a structure in which the names of the fields are given by the
% identifiers in the first column of id_map and the values are given by the
% second column in id_map.  id_map is a cell array.
%
% Example:
%
% id_map = {'EVOKES_EMOTION_ID', 185.1; ...
%      'STRENGTH_EXPER_EMOT_ID', 191.1; ...
%      'EMOTION_TYPE_ID', 192.1};
%
% results in 
%
% s.EVOKES_EMOTION_ID
% s.STRENGTH_EXPER_EMOT_ID
% s.EMOTION_TYPE_ID
%
% where the values in the right column of id_map populate their respective
% fields in s.  There is also a field called all that has all of the IDs in a
% single vector

% 07/31/06 Petr Janata

nid = size(id_map,1);

s.all = [id_map{:,2}];

for iid = 1:nid
  id_str = id_map{iid,1};
  s.(id_str) = id_map{iid,2};
end

return
