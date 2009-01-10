function outdata = ensemble_fmri_mean_hires(indata,defs)

% creates mean hires image from given subjects
% 
% REQUIRES - param.sinfo
% 
% FB 2008.09.10

if (iscell(indata) && ~isempty(indata) && isfield(indata{1},'task') && ...
        ~isempty(strmatch('return_outdir',indata{1}.task))) || ...
        (isstruct(indata) && isfield(indata,'task') && ...
        ~isempty(strmatch('return_outdir',indata.task)))

  outdata = ''; % return nothing, ens_job_|| will use the default
  return
end

global r

r = init_results_struct;

r.type = 'fmri_anal';  % Identify the type of this reporting instance
r.report_on_fly = 1;

outdata = ensemble_init_data_struct();
outdata.type = 'meanhires';

% Parse out the input data
for idata = 1:length(indata)
  if isfield(indata{idata},'type')
    switch indata{idata}.type
      case 'hires'
        hires = indata{idata};
        hicol = set_var_col_const(hires.vars);
    end
  end
end

if isfield(defs,'sinfo') && isstruct(defs.sinfo)
  sinfo = defs.sinfo;
  proc_subs = {sinfo(:).id};
  nsub_proc = length(proc_subs);
end

% check for required vars
check_vars = {'sinfo','hires'};
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

outdata.vars = [outdata.vars 'meanhires'];
mhires_idx = length(outdata.vars);
outdata.data{mhires_idx} = ensemble_init_data_struct();
outdata.data{mhires_idx}.type = 'meanhires';
outdata.data{mhires_idx}.vars = {'path','nsub'};
mhcol = set_var_col_const(outdata.data{mhires_idx}.vars);
outdata.data{mhires_idx}.data{1} = {};
outdata.data{mhires_idx}.data{2} = [];

% get flags
try USE_SPM = defs.mean_hires.USE_SPM; catch USE_SPM = 0; end
try USE_FSL = defs.mean_hires.USE_FSL; catch USE_FSL = 0; end

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
for isub=1:nsub_proc
  subid = sinfo(isub).id;
  msg = sprintf('\tsubject (%d/%d): %s\n', isub, nsub_proc,subid);
  r = update_report(r,msg);
  
  nsess = length(sinfo(isub).sessinfo);
  for isess = 1:nsess
    sess = sinfo(isub).sessinfo(isess);
    
    if ~sess.use_session
	  msg = sprintf('\t\t\tSkipping session %d\n', isess);
	  r = update_report(r,msg);
	  continue
    end
    
    hfilt.include.all.subject_id = {subid};
    hfilt.include.all.session = isess;
    hdata = ensemble_filter(hires,hfilt);
    hpath = hdata.data{hicol.path}{1};
    fprintf('\t\tLoading hires image for session (%d/%d): %s\n',isess,...
        nsess,hpath);
    curr_files = [curr_files; hpath];

  end % for isess =
end % for isub=

if USE_SPM
  % Get file structures
  V = spm_vol(curr_files);

  % Load the data
  Y = spm_read_vols(V);
  Yout = mean(Y,4);

  Vout = V(1);
  Vout.fname = fullfile(defs.paths.groupanat,...
      sprintf('mean_nhires_%dsubs.img',nsub_proc));
  fprintf('Writing group file: %s\n',Vout.fname)

  Vout = spm_write_vol(Vout,Yout);

  outdata.data{mhires_idx}.data{mhcol.path} = [...
      outdata.data{mhires_idx}.data{mhcol.path}; Vout.fname];
  outdata.data{mhires_idx}.data{mhcol.nsub} = [...
      outdata.data{mhires_idx}.data{mhcol.nsub}; nsub_proc];
      
end % if USE_SPM
