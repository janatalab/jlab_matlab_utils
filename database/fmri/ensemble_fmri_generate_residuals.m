function [outdata] = ensemble_fmri_generate_residuals(indata,defs)

% generates residual images for level1 models estimated using SPM
% 
%   outdata = ensemble_fmri_generate_residuals(indata,defs)
% 
% FB 2008.08.27

global defaults r
VERBOSE = 1;

outdata = ensemble_init_data_struct();
outdata.type = 'evaluate contrasts';

r = init_results_struct;

r.type = 'fmri_anal';  % Identify the type of this reporting instance
r.report_on_fly = 1;

% Parse out the input data
for idata = 1:length(indata)
  if isfield(indata{idata},'type')
    switch indata{idata}.type
      case 'sinfo'
        sinfo = indata{idata}.data;
        proc_subs = {sinfo(:).id};
        nsub_proc = length(proc_subs);
      case 'modelspec'
        modelspec = indata{idata};
        mocol = set_var_col_const(modelspec.vars);
      case {'paths'}
        pathdata = indata{idata};
        pcol = set_var_col_const(pathdata.vars);
    end
  end
end

% check for required vars
check_vars = {'sinfo','modelspec','pathdata'};
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
          outdata = spathdata.data{pcol.path}{1};
        else
          sfilt = pfilt;
          sfilt.include.all.path_type = {'sess_outdir'};
          spathdata = ensemble_filter(lpathdata,sfilt);
          if length(spathdata.data{1}) == 1;
            outdata = spathdata.data{pcol.path}{1};
          else
            sfilt = pfilt;
            sfilt.include.all.path_type = {'sub_outdir'};
            spathdata = ensemble_filter(lpathdata,sfilt);
            if length(spathdata.data{1}) == 1;
              outdata = spathdata.data{pcol.path}{1};            
            end
          end
        end
      end
    end
  end
  if ~exist('outdata','var') || ~exist(outdata,'dir'), outdata = ''; end
  return
end

% get flags
try USE_SPM = defs.build_model_l1.USE_SPM; catch USE_SPM = 0; end
try USE_FSL = defs.build_model_l1.USE_FSL; catch USE_FSL = 0; end
try statsdir = defs.build_model_l1.statsdir; catch statsdir = 'stats'; end

if (~USE_FSL && ~USE_SPM) || (USE_FSL && USE_SPM)
  error('you must specify either SPM or FSL to carry out the analyses\n');
end

% FSL not supported yet
if USE_FSL
  error('FSL not yet supported');
end % if USE_SPM

% outdata
% sinfo
outdata.vars = [outdata.vars 'sinfo'];
sinfo_idx = length(outdata.vars);
outdata.data{sinfo_idx} = ensemble_init_data_struct();
outdata.data{sinfo_idx}.type = 'sinfo';
outdata.data{sinfo_idx}.data = sinfo;

% modelspec
outdata.vars = [outdata.vars 'modelspec'];
mod_idx = length(outdata.vars);
outdata.data{mod_idx} = modelspec;

% residual EPI images
outdata.vars = [outdata.vars 'resid_epi'];
res_idx = length(outdata.vars);
outdata.data{res_idx} = ensemble_init_data_struct();
outdata.data{res_idx}.name = 'resid_epi';
outdata.data{res_idx}.type = 'resid_epi';
outdata.data{res_idx}.vars = {'subject_id','session','model_id','path'};
modcol = set_var_col_const(outdata.data{mod_idx}.vars);
outdata.data{res_idx}.data{modcol.subject_id} = {};
outdata.data{res_idx}.data{modcol.session} = [];
outdata.data{res_idx}.data{modcol.model_id} = [];
outdata.data{res_idx}.data{modcol.path} = {};

% 
%  LOOP OVER MODELSPECs
% 
nspec = length(modelspec.data{1});

for j=1:nspec
  model_fname = modelspec.data{mocol.path}{j};
  model_dir = fileparts(model_fname);
  if ~exist(model_fname,'file')
    warning('could not find model spec %s, SKIPPING',model_fname);
    continue
  else
    fprintf(1,'generating residuals for %s\n',model_fname);
  end
  
  % Load the SPM.mat file for the model in order to get the design matrix
  % (X) and list of original files (Y)
  fprintf('Loading model file: %s\n',model_fname);
  load(model_fname);

  % Load the original data (Y)
  fprintf('Loading the original data ...\n');
  Y = spm_read_vols(SPM.xY.VY);

  % Make sure we move to the directory where the beta files live
  start_dir = pwd;
  cd(SPM.swd)

  % Load each of the beta* files (b)
  fprintf('Loading the beta estimates ...\n');
  B = spm_read_vols(SPM.Vbeta);

  % Move back to the directory we started in
  cd(start_dir)

  % Initialize a 4-D output volume, setting it to NaN
  fprintf('Initializing the output data matrix\n');
  res = nan(size(Y));

  % Loop over voxels that are part of the evaluation set
  fprintf('Calculating the residuals ...\n');
  XYZ = SPM.xVol.XYZ;
  nvox = size(XYZ,2);
  for ivox = 1:nvox
    % Extract the vector of beta values for this voxel
    b = squeeze(B(XYZ(1,ivox),XYZ(2,ivox),XYZ(3,ivox),:));

    % Get the residuals for this voxel (res=Y-X*b)
    res(XYZ(1,ivox),XYZ(2,ivox),XYZ(3,ivox),:) = ...
      squeeze(Y(XYZ(1,ivox),XYZ(2,ivox),XYZ(3,ivox),:)) - ...
      SPM.xX.X*b(:);  % the 100 adds an offset to everything so that
    % we don't run into issues downstream
  end % for ivox

  % Add an offset to the residuals such that the minimum value is more than
  % 80% of the mean value
  resmax = nanmax(res(:));
  resmin = nanmin(res(:));
  resmean = nanmean(res(:));
  fprintf('Residuals: Mean=%1.2f, Min=%1.2f, Max=%1.2f\n', resmean, resmin, resmax);

  newmean = max(100,ceil(resmin/(defs.fmri.spm.defaults.mask.thresh-.1-1)));
  fprintf('Adding offset of %d to residuals\n', newmean);
  res = res+newmean;

  %
  % Deal with writing the 4D volume. Parts of this copied from spm_spm.m
  %

  % Make sure we have an output directory for the residuals
  resid_dir = fullfile(model_dir, 'resid'); check_dir(resid_dir);

  % Initialize output images
  fprintf('Initializing output files ...\n');
  nvol = size(res,4);
  clear VResI

  VResI(1:nvol) = deal(struct(...
    'fname',    [],...
    'dim',      SPM.xVol.DIM',...
    'dt',       [spm_type('float64') spm_platform('bigend')],...
    'mat',      SPM.xVol.M,...
    'pinfo',    [1 0 0]',...
    'descrip',  'spm_spm:Residual image'));

  for ivol = 1:nvol
    VResI(ivol).fname   = fullfile(resid_dir,sprintf('ResI_%04d.img', ivol));
    VResI(ivol).descrip = sprintf('spm_spm:ResI (%04d)', ivol);
    spm_unlink(VResI(ivol).fname); % remove any existing copies
  end

  VResI = spm_create_vol(VResI); % initialize the whole image stack

  fprintf('Writing the residuals images ...\n');
  for ivol = 1:nvol
    % write image
    spm_write_vol(VResI(ivol),res(:,:,:,ivol));
    
    % save to outdata struct
    outdata.data{res_idx} = ensemble_add_data_struct_row(...
        outdata.data{res_idx},...
        'subject_id',modelspec.data{mocol.subject_id}{j},...
        'session',   modelspec.data{mocol.session}(j),...
        'model_id',  modelspec.data{mocol.model_id}(j),...
        'path',{VResI(ivol).fname});
  end % for ivol=1:nvol
end % for j=1:nspec
