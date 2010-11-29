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
%       .ncomp - if not specified, icatb_estimate_dimension will be used to
%           estimate ncomp. Note: this might result in a very large # of
%           components NOTE: coded, but not tested
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

global defaults r

%% init data structs
outdata = [];

r = init_results_struct;

r.type = 'fmri_anal';  % Identify the type of this reporting instance
r.report_on_fly = 1;
fprintf(1,'initializing data structures\n');

% convert to cell if necessary
if isstruct(indata)
	indata = {indata};
end

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
if isfield(params.model, 'gift')
	m = params.model.gift;
	m.name = params.model.name;
	m.model_id = params.model.model_id;
else
	warning('PUT YOUR GIFT PARAMETERS IN model.gift DAMMIT!')
	m = params.model;
end

% get number of subjects
subids = unique(epidata.data{epicol.subject_id});
nsub = length(subids);

% Figure out how many sessions we have per subject
numsess_per_sub = nan(nsub,1);
session_ids = cell(nsub,1);
for isub = 1:nsub
	submask = ismember(epidata.data{epicol.subject_id}, subids{isub});
	session_ids{isub} = unique(epidata.data{epicol.session}(submask));
	numsess_per_sub(isub) = length(session_ids{isub});
end
	
nsess = unique(numsess_per_sub);
treatMultipleSessionsAsOne = 0;
if length(nsess) > 1
	fprintf('Found variable numbers of sessions per subjects: %s\n', sprintf('%d ', nsess));
	if isfield(m, 'treatMultipleSessionsAsOne') && m.treatMultipleSessionsAsOne
		treatMultipleSessionsAsOne = 1;
		nsess = 1;
	else
		fprintf('Please resolve this somehow before continuing ...\n');
		fprintf('For example, set treatMultipleSessionsAsOne in model_spec to 1\n');
		return
	end
end

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
try spmMatFlag = m.spmMatFlag;
catch
  if exist('modelspec','var')
    spmMatFlag = 'diff_sub_diff_sess';
  else
    spmMatFlag = 'no';
  end
end
try ncomp = m.ncomp; end
if ~exist('ncomp','var') || isempty(ncomp)
  ncomp=ceil(icatb_estimate_dimension(char(epidata.data{epicol.path})));
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
    numOfGs = [nsub*nsess nsub 1];
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
k = 0;
for isub=1:nsub
  sfilt.include.all.subject_id = subids(isub);
	
	for isess = 1:nsess
		k = k+1;
		
		if isnumeric(session_ids{isub})
			sfilt.include.all.session = session_ids{isub}(isess);
		else
			sfilt.include.all.session = session_ids{isub}{isess};
		end
		sdata = ensemble_filter(epidata,sfilt);
		
		sesInfo.userInput.files(k).name = char(sdata.data{epicol.path});
    sesInfo.userInput.diffTimePoints(k) = length(sdata.data{epicol.path});
  
		% design matrices?
		if exist('modelspec','var')
			mdata = ensemble_filter(modelspec,sfilt);
			sesInfo.userInput.designMatrix(k).name = mdata.data{mcols.path}{1};
		end % if exist('modelspec','var
	end % for isess
end % for k=1:nsub

fprintf(1,'other options\n');

% pca options
if ~exist('pcaOpts','var')
  sesInfo.userInput.pca_opts = icatb_pca_options(pcaType);
else
  sesInfo.userInput.pca_opts = icatb_pca_options(pcaType,pcaOpts,'off');
end

% Get image volume info. Just base this on the first file.
V = spm_vol(epidata.data{epicol.path}{1});
vol_dims = V.dim;

voxdims = diag(V.mat);  % potentially need to deal with issue of first voxdim, x=-1
%sesInfo.userInput.HInfo = struct('DIM',V.dim,'V',V,...
%    'VOX',params.fmri.spm.opts.spmopts.normalise_ropts);  % PJ - isn't this a bit dangerous? Why not pull VOX info from the volume information?
sesInfo.userInput.HInfo = struct('DIM',V.dim,'V',V,'VOX',voxdims(1:3));
	
% Deal with mask file
if isfield(m,'maskFile') && ~isempty(m.maskFile)
	sesInfo.userInput.maskFile = m.maskFile;
	V = spm_vol(sesInfo.userInput.maskFile);
	Y = spm_read_vols(V);
	sesInfo.userInput.mask_ind = find(Y);
else
	fprintf('Constructing EPI mask files');
	
	clear sfilt
	k = 0;
	Ymask_all = nan([vol_dims nsub]);
	for isub = 1:nsub
		sfilt.include.all.subject_id = subids(isub);
	
		Ymask_sub = nan([vol_dims nsess]);
		for isess = 1:nsess		
			if isnumeric(session_ids{isub})
				sfilt.include.all.session = session_ids{isub}(isess);
			else
				sfilt.include.all.session = session_ids{isub}{isess};
			end
			sdata = ensemble_filter(epidata,sfilt);

			% Get file information for this subject/session
			V = spm_vol(char(sdata.data{epicol.path}));
			
			% Load the files
			Y = spm_read_vols(V);
			
			% Calculate the mean
			meanY = mean(Y(:));

			thresh = meanY;  % mean/8 is SPM default
			
			% Create a masked volume
			Ymask = all(Y >= thresh,4); % see which voxels are above threshold for all volumes
			
			% Add it to an overall stack
			Ymask_sub(:,:,:,isess) = Ymask;
			
			% Save the masked volume to a file if desired
			
		end % for isess
		Ymask_all(:,:,:,isub) = all(Ymask_sub,4);
		
	end % for isub
	
	% Get the global list of indices
	sesInfo.userInput.mask_ind = find(all(Ymask_all,4));
	
	% Save the global mask if desired
end

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
