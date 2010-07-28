function outdata = ensemble_fmri_gift_init(indata,params)

% initializes and runs GIFT analyses
% 
%   outdata = ensemble_fmri_gift(indata,params)
% 
% This function will create a GIFT analyis parameter structure, save it to
% disk (in the directory fullfile(params.paths.analroot,'gift',...
% params.model.name)), and run the analysis related to that file.
% 
% NOTE: this currently only allows for one session per participant.
% 
% REQUIRES
%   indata
%       epi_data, realign_epi, resid_epi
%       (modelspec)
%   params.paths.analpath
%   params.fmri.spm.opts.spmopts.normalise_ropts
%   params.model
%       .name
%       .ncomp
%       .algorithm( must specify either a number or a name; see also
%           icatb_icaAlgorithm)
%       a number of parameters may be specified to override defaults for
%       sesInfo.userInput settings:
%           .dataSelMethod (2), .nreduc (2), .preproc_type ('remove mean
%           per timepoint'), .backReconType ('gica3'), .mscaleType (2),
%           .pcaType ('standard'), .pcaOpts (defaults from
%           icatb_pca_options), .spmMatFlag (if no modelspec indata are
%           provided, then this defaults to 'no', otherwise it defaults to
%           'diff_sub_diff_sess')
%       .icasso - if specified (as a struct, with the following
%           fieldnames), ICASSO will be performed
%           .sel_mode, .num_ica_runs
%   params.run_analysis (0 [1] 2:7) - if not set to 0, the given group ICA
%       step (see icatb_runAnalysis.m) will be executed after creating an
%       ica parameter info .mat file from sesInfo
%   
% RETURNS
%   outdata - an ensemble data structure containing the path of the
%       generated sesInfo structure.
% 
% FB 2010.07.06

%% init data structs
outdata = [];

fprintf(1,'initializing data structures\n');

% parse input data
for idata = 1:length(indata)
  if isfield(indata{idata},'type')
      switch indata{idata}.type
        case {'epi','realign_epi','resid_epi'}
          epidata = indata{idata};
          epicol = set_var_col_const(epidata.vars);
        case 'modelspec'
          modelspec = indata{idata};
          mcols = set_var_col_const(modelspec.vars);
      end
  end
end

% check for required vars, quit if they can't be found
check_vars = {'epidata'};
check_required_vars;

if (iscell(indata) && ~isempty(indata) && isfield(indata{1},'task') && ...
        ~isempty(strmatch('return_outdir',indata{1}.task))) || ...
        (isstruct(indata) && isfield(indata,'task') && ...
        ~isempty(strmatch('return_outdir',indata.task)))
  if exist('pathdata','var') && ~isempty(pathdata.data{1})
    if length(nsub_proc) == 1
      pfilt = struct();
      pfilt.include.all.subject_id = proc_subs;
      lpathdata = ensemble_filter(pathdata,pfilt);
      if ~isempty(lpathdata.data{1})
        sfilt = pfilt;
        sfilt.include.all.path_type = {'anal_outdir'};
        spathdata = ensemble_filter(lpathdata,sfilt);
        if length(spathdata.data{1}) == 1
          % one epi outdir, save outdata = epi_outdir
          outdata = spathdata.data{pacol.path}{1};
        else
          sfilt = pfilt;
          sfilt.include.all.path_type = {'sess_outdir'};
          spathdata = ensemble_filter(lpathdata,sfilt);
          if length(spathdata.data{1}) == 1;
            outdata = spathdata.data{pacol.path}{1};
          else
            sfilt = pfilt;
            sfilt.include.all.path_type = {'sub_outdir'};
            spathdata = ensemble_filter(lpathdata,sfilt);
            if length(spathdata.data{1}) == 1;
              outdata = spathdata.data{pacol.path}{1};            
            end
          end
        end
      end
    end
  end
  if ~exist('outdata','var') || ~exist(outdata,'dir'), outdata = ''; end
  return
end

% init outdata
outdata = ensemble_init_data_struct();
outdata.type = 'gift_sesInfo';
if isfield(params,'outDataName')
  outdata.name = params.outDataName;
else
  outdata.name = outdata.type;
end
outdata.vars = {'model_id','path'};
ocol = set_var_col_const(outdata.vars);
outdata.data{ocol.model_id} = [];
outdata.data{ocol.path} = {};

% init sesInfo
sesInfo = struct('userInput',struct(),'isInitialized',0);

%% init vars
fprintf(1,'initializing variables\n');

% model information
m = params.model;

% get number of subjects
subids = unique(epidata.data{epicol.subject_id});
nsub = length(subids);
nsess = 1; % FIXME: parse out sessions (runs) in some way?
nimg = length(epidata.data{1});

% analysis root
ar = fullfile(params.paths.analpath,'gift',m.name);
check_dir(ar,0,1);

% icatb defaults ... FIXME: are defaults appropriate?
try dataSelMethod = m.dataSelMethod; catch dataSelMethod = 2; end
try nreduc = m.nreduc; catch nreduc = 2; end
try preproc_type = m.preproc_type;
catch preproc_type = 'remove mean per timepoint'; end
try pcaType = m.pcaType; catch pcaType = 'standard'; end
try pcaOpts = m.pcaOpts; end
try backReconType = m.backReconType; catch backReconType = 'gica3'; end
try scaleType = m.scaleType; catch scaleType = 2; end
try ncomp = m.ncomp; catch error('ESTIMATE NCOMP'); end
try spmMatFlag = m.spmMatFlag;
catch
  if exist('modelspec','var')
    spmMatFlag = 'diff_sub_diff_sess';
  else
    spmMatFlag = 'no';
  end
end

switch m.algorithm
  case {1,'infomax','InfoMax','infoMax'}
    algo = 1;
  case {2,'FastICA','fastica','fastIca','fastICA'}
    algo = 2;
  otherwise
    error('unknown algorithm\n');
end % switch m.algorithm
icaOpts = icatb_icaOptions([ncomp nimg],algo,'off');

switch nreduc
  case 3
    numOfPC = repmat(ncomp,1,3);
    numOfGs = [nsub*nsess 1 1];
    error('need to think about num of groups 2');
  case 2
    numOfPC = [ncomp ncomp 0];
    numOfGs = [nsub*nsess 1 0];
  case 1
    numOfPC = [ncomp 0 0];
    numOfGs = [nsub*nsess 0 0];
  otherwise
    error('nreduc out of bounds');
end % switch nreduc

% run analysis?
try run_analysis = params.run_analysis; catch run_analysis = 1; end

%% populate sesInfo.userInput
fprintf(1,'populating sesInfo.userInput\n');
sesInfo.userInput = struct(...
    'pwd',ar,...
    'prefix',m.name,...
    'param_file',fullfile(ar,sprintf('%s_ica_parameter_info.mat',m.name)),...
    'dataType','real',...
    'read_complex_images','real&imaginary',...
    'write_complex_images','real&imaginary',...
    'read_complex_file_naming',{{'R_' 'I_'}},...
    'write_complex_file_naming',{{'R_' 'I_'}},...
    'files',struct('name',[]),...
    'dataSelMethod',dataSelMethod,...
    'designMatrix',struct('name',[]),...
    'spmMatFlag',spmMatFlag,...
    'diffTimePoints',zeros(1,nsub*nsess),...
    'modality','fMRI',...
    'numOfSub',nsub,...
    'numOfSess',nsess,...
    'numComp',ncomp,...
    'numReductionSteps',nreduc,...
    'numOfPC1',numOfPC(1),'numOfPC2',numOfPC(2),'numOfPC3',numOfPC(3),...
    'numOfGroups1',numOfGs(1),'numOfGroups2',numOfGs(2),'numOfGroups3',numOfGs(3),...
    'maskFile',[],...
    'preproc_type',preproc_type,...
    'pcaType',pcaType,...
    'pca_opts',struct(),...
    'backReconType',backReconType,...
    'scaleType',scaleType,...
    'HInfo',struct(),...
    'mask_ind',[],...
    'algorithm',algo,...
    'ICA_Options',{icaOpts}...
    );

% files
fprintf(1,'adding files\n');
for k=1:nsub
  sfilt.include.all.subject_id = subids(k);
  sdata = ensemble_filter(epidata,sfilt);
  sesInfo.userInput.files(k).name = cell2mat(sdata.data{epicol.path});
  
  sesInfo.userInput.diffTimePoints(k) = length(sdata.data{epicol.path});
  
  % design matrices?
  if exist('modelspec','var')
    mdata = ensemble_filter(modelspec,sfilt);
    sesInfo.userInput.designMatrix(k).name = mdata.data{mcols.path}{1};
  end % if exist('modelspec','var
end % for k=1:nsub

fprintf(1,'other options\n');

% pca options
if ~exist('pcaOpts','var')
  sesInfo.userInput.pca_opts = icatb_pca_options(pcaType);
else
  sesInfo.userInput.pca_opts = icatb_pca_options(pcaType,pcaOpts,'off');
end

% HInfo - header for first image file, serves as mask
V = spm_vol(epidata.data{epicol.path}{1});
sesInfo.userInput.HInfo = struct('DIM',V.dim,'V',V,...
    'VOX',params.fmri.spm.opts.spmopts.normalise_ropts);

% indices of masked-in voxels
Y = spm_read_vols(V);
sesInfo.userInput.mask_ind = find(Y(:) > mean(Y(:)));

% ICASSO options?
if isfield(m,'icasso')
  sesInfo.userInput.icasso_opts = struct(...
      'sel_mode',m.icasso.sel_mode,'num_ica_runs',m.icasso.num_ica_runs...
      );
end % if isfield(m,'icasso

%% save data
save(sesInfo.userInput.param_file,'sesInfo');
outdata.data{ocol.model_id} = m.model_id;
outdata.data{ocol.path} = {sesInfo.userInput.param_file};

%% execute analysis?
if run_analysis, icatb_runAnalysis(sesInfo,run_analysis); end

%% done!
fprintf(1,'DONE\n');
