function outdata = ensemble_fmri_contrast_mask(indata,defs)

% creates a mask from a given contrast
% 
% function outdata = ensemble_fmri_contrast_mask(indata,defs)
% 
% REQUIRES
% 
% RETURNS
%
% Creates a mask based on a particular contrast.
%
% Currently VERY hard-coded for Tonreg_F
% 
% 10/16/08 PJ Modified to accommodate changes in spm_getSPM.m
% 01/12/09 FB adapted from autobio_contrast_mask.m
% 
% FIXME: currently deals only with one contrast at a time ... could deal
% with creating multiple masks from a given model ...

global r

outdata = ensemble_init_data_struct();
outdata.type = 'contrast_masks';

r = init_results_struct;

r.type = 'fmri_anal';  % Identify the type of this reporting instance
r.report_on_fly = 1;

% Parse out the input data
for idata = 1:length(indata)
  if isfield(indata{idata},'type')
    switch indata{idata}.type
      case {'modelspec','level2model'}
        modelspec = indata{idata};
        mocol = set_var_col_const(modelspec.vars);
    end
  end
end

if isfield(defs,'sinfo') && isstruct(defs.sinfo)
  sinfo = defs.sinfo;
  proc_subs = {sinfo(:).id};
  nsub_proc = length(proc_subs);
end

% check for required vars, quit if they can't be found
check_vars = {'sinfo','modelspec'};
check_required_vars;

% sinfo
outdata.vars = [outdata.vars; 'sinfo'];
sinfo_idx = length(outdata.vars);
outdata.data{sinfo_idx} = ensemble_init_data_struct();
outdata.data{sinfo_idx}.type = 'sinfo';
outdata.data{sinfo_idx}.data = sinfo;

% modelspec
outdata.vars = [outdata.vars 'contrast_masks'];
cm_idx = length(outdata.vars);
outdata.data{cm_idx} = ensemble_init_data_struct();
outdata.data{cm_idx}.type = 'contrast_masks';
outdata.data{cm_idx}.vars = {'model_id','level','subject_id',...
    'contrast','path'};
cmcol = set_var_col_const(outdata.data{cm_idx}.vars);
outdata.data{cm_idx}.data{1} = [];
outdata.data{cm_idx}.data{2} = [];
outdata.data{cm_idx}.data{3} = {};
outdata.data{cm_idx}.data{4} = {};
outdata.data{cm_idx}.data{5} = {};

if (iscell(indata) && ~isempty(indata) && isfield(indata{1},'task') && ...
        ~isempty(strmatch('return_outdir',indata{1}.task))) || ...
        (isstruct(indata) && isfield(indata,'task') && ...
        ~isempty(strmatch('return_outdir',indata.task)))
  if isfield(defs,'paths') && isfield(defs.paths,'analpath')
    outdata = defs.paths.analpath;
    check_dir(outdata);
  end
  if ~exist('outdata','var') || ~exist(outdata,'dir'), outdata = ''; end
  return
end

cp = struct(); % contrast parameters
try cp.thresDesc = defs.contrast_mask.thresDesc;
    catch cp.thresDesc = ''; end % thresDesc = 'FWE';
try cp.u = defs.contrast_mask.thresh; catch cp.u = 0.005; end
try cp.k = defs.contrast_mask.extent; catch cp.k = 1; end
try cp.Im = defs.contrast_mask.mask; catch cp.Im = []; end
try cp.title = defs.contrast_mask.title; catch cp.title = ''; end
try LEVEL1 = defs.contrast_mask.LEVEL1; catch LEVEL1 = 0; end
try LEVEL2 = defs.contrast_mask.LEVEL2; catch LEVEL2 = 0; end

msg = '';
if isfield(defs,'contrasts') && ~isempty(defs.contrasts) && ...
        iscell(defs.contrasts) && ~isempty(defs.contrasts)
  contrasts = defs.contrasts;
  ncont = length(contrasts);
else
  msg = 'please specify a cell array of contrasts in defs.contrasts\n';
end

if isfield(defs,'model') && ~isempty(defs.model)
  model = defs.model;
  model_id = model.model_id;
else
  msg = 'please specify a contrast in defs.cont_name\n';
end

if ~LEVEL1 && ~LEVEL2
  msg = [msg 'must provide a model level or levels for contrasts\n'];
end

if LEVEL1 && nsub_proc < 1
  msg = [msg 'must provide subject data if making level1 contrast masks\n'];
end

if ~isempty(msg), error(msg), end

for icont=1:ncont
  cp.cont_name = contrasts{icont};

  if LEVEL1
    sf.include.model_id = model_id;
  
    for isub = 1:nsub_proc
      sf.include.subject_id = proc_subs(isub);
    
      % Determine number of sessions for this subject
      nsess = length(sinfo(isub).sessinfo);
  
      %
      % START OF THE SESSION LOOP
      %
  
      for isess = 1:nsess

        ssf = sf;
        ssf.include.session = isess;
        ssmod = ensemble_filter(modelspec,ssf);
        model_fname = ssmod.data{mocol.path}{1};
      
        %
        % Set up the contrast that we will use to generate the mask
        %
        load(model_fname);
        cp.xCon = SPM.xCon;
        cp.swd = fileparts(model_fname); % sanity check?
        cp.Ic = strmatch(cont_name, {SPM.xCon.name}, 'exact'); % sanity check?
      
        cname = eval_contrast(cp);

        outdata.data{cm_idx}.data{cmcol.model_id} = [...
            outdata.data{cm_idx}.data{cmcol.model_id}; model_id];
        outdata.data{cm_idx}.data{cmcol.level} = [...
            outdata.data{cm_idx}.data{cmcol.level}; 1];
        outdata.data{cm_idx}.data{cmcol.subject_id} = [...
            outdata.data{cm_idx}.data{cmcol.subject_id}; proc_subs{isub}];
        outdata.data{cm_idx}.data{cmcol.contrast} = [...
            outdata.data{cm_idx}.data{cmcol.contrast}; cp.cont_name];
        outdata.data{cm_idx}.data{cmcol.path} = [...
            outdata.data{cm_idx}.data{cmcol.path}; cname];
      
      end % for isess = 1:nsess
    end % for isub = 1:nsub_proc
  end % if LEVEL1

  if LEVEL2
    sf.include.model_id = model_id;
    sf.include.contrast = {cont_name};
    smod = ensemble_filter(modelspec,sf);
    model_fname = smod.data{mocol.path}{1};

    %
    % Set up the contrast that we will use to generate the mask
    %
    load(model_fname);
    cp.xCon = SPM.xCon;
    cp.swd = fileparts(model_fname); % sanity check?
    cp.Ic = strmatch(cont_name, {SPM.xCon.name}, 'exact'); % sanity check?

    cname = eval_contrast(cp);
    
    outdata.data{cm_idx}.data{cmcol.model_id} = [...
        outdata.data{cm_idx}.data{cmcol.model_id}; model_id];
    outdata.data{cm_idx}.data{cmcol.level} = [...
        outdata.data{cm_idx}.data{cmcol.level}; 1];
    outdata.data{cm_idx}.data{cmcol.subject_id} = [...
        outdata.data{cm_idx}.data{cmcol.subject_id}; proc_subs{isub}];
    outdata.data{cm_idx}.data{cmcol.contrast} = [...
        outdata.data{cm_idx}.data{cmcol.contrast}; cp.cont_name];
    outdata.data{cm_idx}.data{cmcol.path} = [...
        outdata.data{cm_idx}.data{cmcol.path}; cname];

  end % if LEVEL2
end % for icont=1:ncont

% % % % 
% % % SubRoutines
% % % % 

% function eval_contrast(cp)

function cname = eval_contrast(cp)

% Evaluate the contrast
[SPM,xSPM] = spm_getSPM(cp);
nvox = size(xSPM.XYZ,2);

% Load the current mask volume just to get an SPM volume template
Vorig = spm_vol(fullfile(cp.swd,'mask.img'));
Vorig = rmfield(Vorig,'private');

% Create the mask volume
Vmask = Vorig;
Ymask = zeros(Vmask.dim);
for ivox = 1:nvox
  Ymask(xSPM.XYZ(1,ivox),xSPM.XYZ(2,ivox),xSPM.XYZ(3,ivox)) = 1;
end

cname = fullfile(cp.swd,sprintf('%s_mask.img', cp.cont_name));
if exist(cname,'file'), delete(cname), end
Vmask.fname = cname;
Vmask.descrip = sprintf('%s: threshtype=%s, thresh=%1.4f, extent=%d', ...
    cp.cont_name,cp.thresDesc, cp.u, cp.k);
fprintf('Writing contrast mask: %s\n', Vmask.fname);
Vmask = spm_write_vol(Vmask,Ymask);
