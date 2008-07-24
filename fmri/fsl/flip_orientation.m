function flip_orientation(infname, outfname)
% flip_orientation(infname)
%
% Flips the input image from radiological to neurological or vice versa
%
% 

% 12/02/05 Petr Janata

% Swap the data
fsl_str = sprintf('/usr/local/fsl/bin/avwswapdim %s -x y z %s', infname, outfname);
fprintf('%s\n', fsl_str);
status = unix(fsl_str);
if status
  error(sprintf('Failed to swap data orientation <%s>\n', infname))
end
    
% Swap the header
fsl_str = sprintf('/usr/local/fsl/bin/avworient -swaporient %s', outfname);
fprintf('%s\n', fsl_str);
status = unix(fsl_str);
if status
  error(sprintf('Failed to swap data orientation <%s>\n', outfname))
end
