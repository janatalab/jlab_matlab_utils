function [names,vals] = fmri_regress_miremotion(pinfo,minfo,sess)

% generates regressors from MIRToolbox miremotion output
% 
%   [names,vals] = fmri_regress_miremotion(pinfo,minfo,sess)
% 
% this script retrieves data from previously-calculated mirtoolbox function
% miremotion, uses spline interpolation to resample the data into fmri-ready
% temporal space (TR/dt), and then returns a matrix of regressors.
% 
% REQUIRES
%   pinfo.mysql (containing .host, .database, .user, .passwd, .conn_id
%   pinfo.scanner.(dt/TR/actual_nvol)
%   minfo.response_filter
%   minfo.music_dur_max
%   minfo.mirtoolbox.miremotion
%       .params - parameters that match the calculated features that you
%           would like to find on disk
%       .stime - length of each sample in seconds
%       .mirgetdata - set to 1 if you want to take the mir function output
%           and extract regressor signal using mirgetdata()
%       .get_data - cell array of strings identifying get_data arguments
%           that will be used to extract data from a given function output.
%           If mirgetdata is not used here, specifying get_data will direct
%           this script to the data that you want to extract.
%       .get_proc_fh - if specified, use_data will be filtered through the
%           function specified in this field before being further
%           processed.
%   minfo.leaky_integrator - if specified as a struct with the following
%       fields, data from mirtoolbox will be passed through
%       IPEMLeakyIntegration before being interpolated
%       .fs - sampling frequency of the given input data
%       .HalfDecayTime - time (in s) at which an impulse would be reduced
%           to half its value
% 
%   mirgetdata supercedes get/use_data, but at least one must be specified
% 
%   if DimData or ClassData are specified in get_data, this function will
%   output separate regressors for each of the dims or classes, and will
%   specify the regressor names as 'miremotion%Dim%' or 'miremotion%Class'
%   (i.e. miremotionActivity ... miremotionAnger)
% 
% FB 2010.06.27

%% init vars
names = {};
vals = [];

ensemble_globals;

m = pinfo.mysql;

TR = pinfo.scanner.TR;
dt = pinfo.scanner.dt;
actual_nvol = pinfo.scanner.actual_nvol;
npts = actual_nvol*dt;

vals = zeros(actual_nvol,0);

% get mirtoolbox info
mirinfo = minfo.mirtoolbox;
mirfuncs = fieldnames(mirinfo);
nmirfunc = length(mirfuncs);

% Get a list of the IDs to process
sfilt.include.all = minfo.response_filter;
sdata = ensemble_filter(pinfo,sfilt);
pc = set_var_col_const(pinfo.vars);

stim_ids = cellfun(@str2num,sdata.data{pc.EVENT_CODE});
nids = length(stim_ids);
onsets = sdata.data{pc.RUN_REL_TIME}/1000;

if ~nids
  fprintf(1,'no sound stimuli in this run, SKIPPING\n');
  return
end

% init funcname struct
fnstr.DimData = {'Activity','Valence','Tension'};
fnstr.ClassData = {'Happy','Sad','Tender','Fear','Anger'};
fnstr.ActivityFactors = {'ActF1','ActF2','ActF3','ActF4','ActF5'};
fnstr.ValenceFactors  = {'ValF1','ValF2','ValF3','ValF4','ValF5'};
fnstr.TensionFactors  = {'TnsF1','TnsF2','TnsF3','TnsF4','TnsF5'};
fnstr.HappyFactors    = {'HapF1','HapF2','HapF3','HapF4','HapF5'};
fnstr.SadFactors      = {'SadF1','SadF2','SadF3','SadF4','SadF5'};
fnstr.TenderFactors   = {'TndF1','TndF2','TndF3','TndF4','TndF5'};
fnstr.FearFactors     = {'FeaF1','FeaF2','FeaF3','FeaF4','FeaF5'};
fnstr.AngerFactors    = {'AngF1','AngF2','AngF3','AngF4','AngF5'};

% get stim names and paths
fprintf(1,'getting stimulus paths\n');

stim_paths = cell(nids,2); % (:,1) - name, (:,2) - path
for k=1:nids
  stim_id = stim_ids(k);

  % get stimulus path
  floc = mysql_extract_data('conn_id',m.conn_id,...
      'table','stimulus', ...
      'extract_flds','location', ...
      'stimulus_id',stim_id);
  [stim_paths{k},stim_names{k}] = fileparts(floc{1}{1});
end

% % % 
% % %   we should handle this, but not at the moment
% % % 
% % are we segregating regressors?
% if isfield(minfo,'tonreg')
%   % segregation params
%   tregp = minfo.tonreg;
%   segnames = tregp.seg_names;
%   segvals = tregp.seg_vals;
%   nseg = length(segvals);
%   ccm = parse_fh(minfo.cond_cue_map);
%   cue = ccm(tregp.seg_cue);
%   resp_params = extract_resp_params_v2(cue,pinfo,minfo,sess);
% else
%   tregp = '';
% end

% get data from each function for all stims
fprintf(1,'Generating regressors\n');
for k=1:nmirfunc
  fname = mirfuncs{k};
  if ~isfield(mirinfo.(fname),'params')
    fparams = struct();
  else
    fparams = mirinfo.(fname).params;
  end
  
  % preproc function?
  if isfield(mirinfo.(fname),'get_proc_fh')
    get_proc_fh = parse_fh(mirinfo.(fname).get_proc_fh);
  end
  
  % initialize data matrix
  if isfield(mirinfo.(fname),'get_data')
    get_data = mirinfo.(fname).get_data;
    funcnames = {};
    for l=1:length(get_data)
      funcnames = [funcnames fnstr.(get_data{l})];
    end % for l=1:length(get_data)
    funcnames = strcat(fname,funcnames);
    nget = length(funcnames);
    funcsig = zeros(npts,nget);
  else
    funcsig = zeros(npts,1);
    funcnames = {fname};
  end
  
  % iterate over stims
  for l=1:nids
    stim_id = stim_ids(l);
    fprintf(1,'feature %s, stim id %d\n',fname,stim_id);

    % search for matching file
    fpath = fullfile(stimulus_root,stim_paths{l},fname);
    fdirs = dir(fullfile(fpath,sprintf('*%s*mat',fname)));
    
    match = 0;
    for m=1:length(fdirs)
      if fdirs(m).isdir, continue, end
      clear job_data;
      load(fullfile(fpath,fdirs(m).name));
      if ~exist('job_data','var'), continue, end
      if ~isfield(job_data,'params') || ~isstruct(job_data.params) || ...
              ~isfield(job_data.params,'mirtoolbox') || ...
              ~isstruct(job_data.params.mirtoolbox) || ...
              ~isfield(job_data.params.mirtoolbox,'params') || ...
              ~isstruct(job_data.params.mirtoolbox.params)
        job_data.params.mirtoolbox.params = struct();
      end
      if compare_structs(job_data.params.mirtoolbox.params,fparams) && ...
              compare_structs(fparams,job_data.params.mirtoolbox.params)
        if isfield(job_data,'results') && ~isempty(job_data.results)
          match = 1;
          break
        end
      end
    end % for m=1:length(fdirs
    
    if ~match
%       error('could not find %s for stimulus id %d',fname,stim_id);
        warning('could not find %s for stimulus id %d, SKIPPING',fname,stim_id);
        continue
    end
    
    % get data
    if ~isempty(strmatch(job_data.results.type,fname,'exact'))
      fdata = job_data.results;
    else
      idx = strmatch(fname,job_data.results.vars,'exact');
      if isempty(idx)
        error('could not find %s in results from %s for stimulus id %d',...
            fname,fdirs(m).file,stim_id);
      end
      fdata = job_data.results.data{idx(1)};
    end
    fcol = set_var_col_const(fdata.vars);

    fsig = [];
    
    % initialize the creator function object
    try
      a = miraudio(rand(1,10000));
      e = miremotion(a);
    end
      
    % load data
    mirdata = load(fdata.data{fcol.mirout_path}{1});
      
    if isfield(mirinfo.(fname),'mirgetdata')
      % use mirgetdata to get data
      fsig = mirgetdata(mirdata.a);
    elseif isfield(mirinfo.(fname),'get_data')
      for m=1:length(get_data)
        lget = get_data{m};
        ldata = get(mirdata.a,lget);
        if isfield(mirinfo.(fname),'get_proc_fh')
          ldata = get_proc_fh(ldata{1});
        end
        fsig = [fsig ldata];
      end % for m=1:nuse
    end % if isfield(mirinfo.(fname
    
    % leaky integrator?
    if isfield(minfo,'leaky_integrator') && ~isempty(minfo.leaky_integrator)
      leakfs = minfo.leaky_integrator.fs;
      leakdt = minfo.leaky_integrator.HalfDecayTime;
      fsig = IPEMLeakyIntegration(fsig',leakfs,leakdt)';
    end % if isfield(minfo,'leaky_integrator

    % spline interpolation!
    stime = mirinfo.(fname).stime; %%% need to make sure to specify in model def!!
    slen  = size(fsig,1)*stime;
    svals = 0:stime:(slen - stime);
    ovals = 0:TR/dt:(slen - TR/dt);
    ssig  = [];
    for m=1:size(fsig,2)
      ssig = [ssig spline(svals,fsig(:,m),ovals)'];
    end % for m=1:size(fsig
    
    % cut off at music_dur_max?
    if isfield(minfo,'music_dur_max')
       % make sure that curr_val doesn't exceed music_dur_max
      dmax = minfo.music_dur_max*dt/TR;
      if size(ssig,1) > dmax
        ssig(dmax+1:end,:) = [];
      end
    end

    onset = round(onsets(l)*dt/TR);
    ofset = onset + size(ssig,1) - 1;
    funcsig(onset:ofset,:) = ssig;
  end % for l=1:nids
  
  % convolve with HRF, add to vals
  convreg = [];
  for l=1:size(funcsig,2)
    convreg(:,l) = conv(funcsig(:,l),spm_hrf(TR/dt));
    vals = [vals convreg([0:actual_nvol-1]*dt+1,l)];
  end % for l=1:size(funcsig,2
  
  names = [names funcnames];
end % for k=1:nmirfunc

% scale to ~ -1:1
% vals = vals - mean(vals(:));
% vals = vals./max(max(vals(:)),abs(min(vals(:))));
% vals = vals - repmat(mean(vals),size(vals,1),1);
