function st = convert_charfld2cellstrfld(st)
% Converts character array fields of a structure to cell arrays of strings
%
% st = convert_charfld2cellstrfld(st);

% 18Nov2012 Petr Janata

fldnames = fieldnames(st);
nfld = length(fldnames);
for ifld = 1:nfld
  currFld = fldnames{ifld};
  if ischar(st.(currFld))
    st.(currFld) = cellstr(st.(currFld));
  end
end

end