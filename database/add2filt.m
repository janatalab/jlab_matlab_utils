function outfilt = add2filt(oldfilt,newfilt)
% Adds filtering fields in the newfilt structure to the oldfilt structure.
%
% outfilt = add2filt(oldfilt,newfilt);
%
% These are the filtering fields that are used by ensemble_apply_crit()
%

% 01/26/07 Petr Janata
% 02/21/07 PJ - turned into a recursive function
% 05/02/07 PJ - fixed a bug that prevented concatenation of old and new values
%               for the same fieldname. Thanks to Jason Golubock.
%          JG - add validation check for empty input structs

if( ~length( newfilt ) )
	outfilt = oldfilt;
	return;
end

if( ~length( oldfilt ) )
	outfilt = newfilt;
	return;
end

outfilt = oldfilt;

old_flds = fieldnames(oldfilt);
new_flds = fieldnames(newfilt);

for inew = 1:length(new_flds)
  currfldname = new_flds{inew};
  oldidx = strmatch(currfldname,old_flds);
  if isempty(oldidx)
    outfilt.(currfldname) = newfilt.(currfldname);
  else
    % Check to see if we are dealing with a struct. If so, compare the
    % sub-fields
    if isstruct(newfilt.(currfldname)) && isstruct(oldfilt.(currfldname))
      outfilt.(currfldname) = add2filt(oldfilt.(currfldname),newfilt.(currfldname));
    elseif ~isstruct(newfilt.(currfldname)) & ~isstruct(oldfilt.(currfldname))
      % Combine field values
      new_vals = ...
	  newfilt.(currfldname)(~ismember(newfilt.(currfldname),oldfilt.(currfldname)));
      old_vals = oldfilt.(currfldname);
      outfilt.(currfldname) = [old_vals new_vals];
    else
      fprintf('%s: Field type mismatch in field: %s\n',  mfilename, currfldname);
    end
  end % if isempty(oldidx)
end % for inew
