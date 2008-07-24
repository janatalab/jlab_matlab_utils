function ev_info = make_motion_evs(par_fname,outdir,outstub)
% Parses the 6 movement parameters in a .par file into separate files so that
% each can be loaded as an explanatory variable (EV) in an FSL model.
%
% par_fname - name of the parameter file
% outdir - name of directory into which to place the EV files
% ev_info - a structure containing information about the 6 EVs. This structure
%           is used by other, e.g. model building, scripts

% 03/20/06 Petr Janata

if nargin < 2
  outdir = '';
  outstub = '';
end

ev_info = create_ev_info;

if ~exist(par_fname)
  fprintf('Movement parameter file <%s> does not exist\n', par_fname);
  return
end

% Load the file
data = load(par_fname);

ncol = size(data,2);
if ncol ~= 6
  fprintf('Parameter file <%s> does not have 6 columns\n', par_fname);
  return
end

ev_info.set_type = 'motion';
ev_info.ev_names = {'pitch','roll','yaw','trans_x','trans_y','trans_z'};
ev_info.nev = 6;

for icol = 1:ncol
  outfname = fullfile(outdir,sprintf('%s_motion_ev%d.txt', outstub,icol));
  tmp = data(:,icol);
  
  fprintf('Writing %s\n', outfname);
  save(outfname,'tmp','-ascii');  % save the actual file
  
  % copy the name of the file to the fnames field
  ev_info.fnames{icol} = outfname;
end

return

