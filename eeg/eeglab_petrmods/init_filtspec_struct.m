function fs = init_filtspec_struct
% fs = init_filtspec_struct;
%
% Initializes fields in a structure with filtering parameters
%

fs = struct('low_cutoff',[], ...  % low frequency cutoff
    'high_cutoff', [], ...
    'filt_order', [], ...
    'filt_type', '', ...
    'inpath', '', ...  % directory with input files
    'outpath', '' ...  % directory for output files
    );

return
