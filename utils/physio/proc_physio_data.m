function [physio] = proc_physio_data(flist,p)
% [physio] = proc_physio_data(flist,params,physio);
%
% See also init_physmon_params.m for specific parameters and default values
% associated with a particular physiological data acquisition system
%

% 07/20/06 Petr Janata

if ~isfield(p,'source')
  error('The source of the data must be specified\n');
end

switch p.source
  case 'keithley'
    physio = read_keithley(flist,p);
    
  case 'mate'
    physio = read_mate(flist,p);
    
  case 'biopac'
    physio = read_biopac(flist,p);
    
  otherwise
    fprintf('Unknown source of physio data: %s\n', p.source);
end
