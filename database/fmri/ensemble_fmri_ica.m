function outdata = ensemble_fmri_ica(indata,defs)

% performs individual independent component analyses on given 4D volumes
%
% CURRENTLY only uses FSL's Melodic
% CURRENTLY expects 4D nifti volumes
% CURRENTLY runs separate ICA on each 4D volume provided ...
% 
% REQUIRES
%   epi data - 4D nifti files - it will iterate over each 4D nifti file
%       given, and perform ICA analyses on each file
%   path data
% 
%   defs.filt - filter criteria applied to epi data
%   defs.ica.NDIMS - number of IC dimensions to retain, 0 = default melodic
%       estimation
%   defs.ica.BGTHRESH
%   defs.ica.MMTHRESH
%   defs.ica.TR
%   defs.ica.melodic_params - extra parameters t pass to melodic
%       default: --nobet --Ostats
% 
% output directories:
%   if a given 4D volume has
% 
% 2009.03.05 FB

global defaults r
VERBOSE = 1;

outdata = ensemble_init_data_struct();
outdata.type = 'ica';

r = init_results_struct;

r.type = 'fmri_ica';  % Identify the type of this reporting instance
r.report_on_fly = 1;

dstr = datestr(now(),30);

% Parse out the input data
for idata = 1:length(indata)
  if isfield(indata{idata},'type')
    switch indata{idata}.type
      case {'epi','realign_epi'}
        epidata = indata{idata};
        epicol = set_var_col_const(epidata.vars);
      case {'paths'}
        pathdata = indata{idata};
        pcol = set_var_col_const(pathdata.vars);
    end
  end
end

% check for required vars
check_vars = {'pathdata','epidata'};
check_required_vars;

% return an output directory for ensemble_jobman_parallel
% for this function, it should be the analysis output directory
if (iscell(indata) && ~isempty(indata) && isfield(indata{1},'task') && ...
        ~isempty(strmatch('return_outdir',indata{1}.task))) || ...
        (isstruct(indata) && isfield(indata,'task') && ...
        ~isempty(strmatch('return_outdir',indata.task)))
  if isfield(defs,'paths') && isfield(defs.paths,'analpath')
    outdata = defs.paths.analpath;
  end
  if ~exist('outdata','var') || ~exist(outdata,'dir'), outdata = ''; end
  return
end

% outdata
outdata.vars = [outdata.vars 'ica_reports'];
ica_idx = length(outdata.vars);
outdata.data{ica_idx} = ensemble_init_data_struct();
outdata.data{ica_idx}.type = 'ica_report';
outdata.data{ica_idx}.vars = {'subject_id','session','ensemble_id',...
    'run','path'};
outdata.data{ica_idx}.data{1} = {};
outdata.data{ica_idx}.data{2} = [];
outdata.data{ica_idx}.data{3} = [];
outdata.data{ica_idx}.data{4} = [];
outdata.data{ica_idx}.data{5} = {};

% get flags
try USE_SPM = defs.ica.USE_SPM; catch USE_SPM = 0; end
try USE_FSL = defs.ica.USE_FSL; catch USE_FSL = 1; end
try MMTHRESH = defs.ica.MMTHRESH; catch MMTHRESH = 1; end
try NDIMS = defs.ica.NDIMS; catch NDIMS = 0; end
try TR = defs.ica.TR; end
if ~exist('TR','var')
  try TR = defs.fmri.protocol.epi.tr; catch error('must specify a TR'); end
end
try BGTHRESH = defs.ica.BGTHRESH; catch BGTHRESH = 1; end
try melodic_params = defs.ica.melodic_params;
    catch melodic_params = '--nobet --Ostats';
    end

if USE_SPM && ~USE_FSL
  error('SPM not supported yet ...\n');
elseif ~USE_FSL && ~USE_SPM
  error(['\t\tyou must specify either SPM or FSL to carry out '...
      'the analyses\n']);
end

%
% apply filtering criteria
%

if isfield(defs,'filt')
  epidata = ensemble_filter(epidata,defs.filt);
end

%
% LOOP OVER IMAGES, estimate
%

nimg = length(epidata.data{epicol.path});
opfilt.include.all.path_type = {'anal_outdir'};

for ii=1:nimg

  % get path to image for estimation
  imgpath = epidata.data{epicol.path}{ii};
  [fp,fn,fx] = fileparts(imgpath);
  
  % get output path - try anal_outdir for subject/session from pathdata
  opfilt.include.all.subject_id = epidata.data{epicol.subject_id}(ii);
  opfilt.include.all.session    = epidata.data{epicol.session}(ii);
  opdata = ensemble_filter(pathdata,opfilt);
  
  if ~isempty(opdata.data{pcol.path}) && exist(opdata.data{pcol.path}{1},'dir')
    % 'anal_outdir/ica' for subject/session
    opath = fullfile(opdata.data{pcol.path}{1},'ica');
    check_dir(opath);
    opath = fullfile(opath,fn);
    check_dir(opath);
    opath = fullfile(opath,dstr);
    check_dir(opath);
  else
    % construct from defs.paths.outroot
    opath = defs.paths.outroot;
    check_dir(opath);
    if ~isempty(epidata.data{epicol.subject_id}{ii})
      % outroot/subject_id
      opath = fullfile(opath,epidata.data{epicol.subject_id}{ii});
      check_dir(opath);
      if ~isempty(epidata.data{epicol.session}(ii))
        % outroot/subject_id/sessionN
        opath = fullfile(opath,...
            sprintf('session%d',epidata.data{epicol.session}(ii)));
        check_dir(opath);
      end
      % add 'analyses'
      opath = fullfile('analyses');
      check_dir(opath);
    else
      % outroot/analyses
      opath = defs.paths.analpath;
      check_dir(opath);
    end
    % add 'ica/dstr'
    opath = fullfile(opath,'ica',dstr);
    check_dir(opath);
  end
  
  % run ICA analysis
  if exist(imgpath,'file') || exist([imgpath '.nii'],'file') || ...
          exist([imgpath '.nii.gz'],'file')
    mstr = sprintf(['melodic -i %s -o %s --bgthreshold=%d --tr=%1.1f '...
        '-d %d --mmthresh=%1.1f --report -v %s'],...
        imgpath,opath,BGTHRESH,TR,NDIMS,MMTHRESH,melodic_params);
    status = unix(mstr);
    if status
      warning('melodic failed to analyze %s:\n\t(%s)\n',imgpath,mstr);
    end
  else
    warning('%s not found, no ICA estimated\n',imgpath);
  end
  
end
