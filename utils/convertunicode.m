function outstr = convertunicode(instr)
% Hack conversions of unicode characters

% 27Jun2014 Petr Janata
outstr = instr;
done = 0;

while ~done
  if any(unicode2native(outstr)==239)
    outstr(unicode2native(outstr)==239) = 101;
  else
    done = 1;
  end
end

return
