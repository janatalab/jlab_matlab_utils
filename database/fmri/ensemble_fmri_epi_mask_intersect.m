function outdata = ensemble_fmri_epi_mask_intersect(indata,defs)

% builds EPI intersect mask
% 
% APPARENTLY ASSUMES NORMED EPI MASK? maybe ease it up a little bit to
% accept whatever mask it's given? maybe indata(idata).type should not
% necessarily be "normed mask" but maybe "epi_mask" or
% {'normed_epi_mask','mean_epi_mask','epi_mask'} ... ?
% 
% FB 2008.08.27

if (iscell(indata) && ~isempty(indata) && isfield(indata{1},'task') && ...
        ~isempty(strmatch('return_outdir',indata{1}.task))) || ...
        (isstruct(indata) && isfield(indata,'task') && ...
        ~isempty(strmatch('return_outdir',indata.task)))

  outdata = ''; % return nothing, ens_job_|| will use the default
  return
end

global r;

outdata = ensemble_init_data_struct();
outdata.type = 'epi_mask_intersect';

r = init_results_struct;

r.type = 'fmri_anal';  % Identify the type of this reporting instance
r.report_on_fly = 1;

% Parse out the input data
for idata = 1:length(indata)
  if isfield(indata{idata},'type')
    switch indata{idata}.type
      case 'epi_mask'
        emdata = indata{idata};
        emcol = set_var_col_const(emdata.vars);
    end
  end
end

if isfield(defs,'sinfo') && isstruct(defs.sinfo)
  sinfo = defs.sinfo;
end

% check for required vars
check_vars = {'sinfo','emdata'};
check_required_vars;

if (iscell(indata) && ~isempty(indata) && isfield(indata{1},'task') && ...
        ~isempty(strmatch('return_outdir',indata{1}.task))) || ...
        (isstruct(indata) && isfield(indata,'task') && ...
        ~isempty(strmatch('return_outdir',indata.task)))
  if isfield(defs,'paths') && isfield(defs.paths,'groupanat')
    outdata = defs.paths.groupanat;
    check_dir(outdata);
  end
  if ~exist('outdata','var') || ~exist(outdata,'dir'), outdata = ''; end
  return
end

% outdata
outdata.vars = [outdata.vars 'sinfo'];
sinfo_idx = length(outdata.vars);
outdata.data{sinfo_idx} = ensemble_init_data_struct();
outdata.data{sinfo_idx}.type = 'sinfo';
outdata.data{sinfo_idx}.data = sinfo;

outdata.vars = [outdata.vars 'epi_mask_intersect'];
emi_idx = length(outdata.vars);
outdata.data{emi_idx} = ensemble_init_data_struct();
outdata.data{emi_idx}.type = 'epi_mask_intersect';
outdata.data{emi_idx}.vars = {'path'};
outdata.data{emi_idx}.data{1} = {};

try USE_SPM = defs.epi_mask_intersect.USE_SPM; catch USE_SPM = 0; end
try USE_FSL = defs.epi_mask_intersect.USE_FSL; catch USE_FSL = 0; end

if USE_FSL && ~USE_SPM
  msg = sprintf('FSL not supported yet ...\n');
  r = update_report(r,msg);
  return
elseif ~USE_FSL && ~USE_SPM
  msg = sprintf(['\t\tyou must specify either SPM or FSL to carry out '...
      'the analyses\n']);
  r = update_report(r,msg);
  return
end

curr_files = [];
VS = {};

nsub_proc = length({sinfo(:).id});
for isub=1:nsub_proc
  subid = sinfo(isub).id;

  if strmatch(subid,emdata.data{emcol.subject_id})
    fprintf('\t\tLoading EPI mask image for (%d/%d): %s\n', isub, nsub_proc,subid);

    emFilt.include.all.subject_id = {subid};
    mpath = ensemble_filter(emdata,emFilt);
    curr_files = [curr_files; mpath.data{emcol.path}{1}];
  else
    fprintf('\t\tEPI mask image for %s not found, skipping ...\n\n',subid);
  end

end % for isub=

% Get file structures
V = spm_vol(strvcat(curr_files));

% Load the data
Y = spm_read_vols(V);

% Perform the desired calculation on the data
Yout = all(Y,4);

% Create an output volume
Vout = V(1);
outfile = fullfile(defs.paths.groupanat,sprintf('epi_mask_intersect.img'));

Vout.fname = outfile;
fprintf('Writing group file: %s\n',Vout.fname)
  
Vout = spm_write_vol(Vout,Yout);

outdata.data{emi_idx}.data{1} = {outfile};
