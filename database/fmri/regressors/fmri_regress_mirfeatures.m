function [names,vals] = fmri_regress_mirfeatures(pinfo,minfo,sess)

% generates regressors from MIRToolbox features
% 
%   [names,vals] = fmri_regress_mirfeatures(pinfo,minfo,sess)
% 
% this script retrieves data from previously-calculated mirtoolbox
% functions, uses spline interpolation to resample the data into fmri-ready
% temporal space (TR/dt), and then returns a matrix of regressors.
% 
% REQUIRES
%   pinfo.mysql (containing .host, .database, .user, .passwd, .conn_id
%   pinfo.scanner.(dt/TR/actual_nvol)
%   minfo.response_filter
%   minfo.music_dur_max
%   minfo.mirtoolbox
%       .(feature function).params - parameters that match the calculated
%           features that you would like to find on disk
%       .(feature function).stime - length of each sample in seconds
%       .(feature function).mirgetdata - set to 1 if you want to take the
%           mir function output and extract regressor signal using
%           mirgetdata()
%       .(feature function).get_data - cell array of strings identifying
%           get_data arguments that will be used to extract data from a
%           given function output. If mirgetdata is not used here,
%           specifying get_data will direct this script to the data that
%           you want to extract.
%       .(feature function).get_proc_fh - if specified, use_data will be
%           filtered through the function specified in this field before
%           being further processed.
% 
%       .(feature function).reg_types - cell array of strings indicating
%           what types of regressors to return. legal values include:
%           'mean' - the mean amplitude of the given feature across each
%               song will be modeled with a boxcar regressor with all
%               positive values
%           'mean_zero' - the above regressor (for 'mean') will be created,
%               and zero-centered across songs.
%       	'zvar_within_song' - the signal of each feature will be
%               standardized within song
%       	'zvar_across_song' - the signal of each feature will be
%               standardized across songs
%           'raw' - the raw signal for the feature will be used to generate
%               a regressor.
% 
%           NOTE: all regressors (except for boxcar) will be resampled into
%           the temporal space used by spm for analyses (TR/dt) using
%           spline interpolation, and then convolved with the canonical hrf
% 
%   minfo.leaky_integrator - if specified as a struct with the following
%       fields, data from mirtoolbox will be passed through
%       IPEMLeakyIntegration before being interpolated. NOTE: leak
%       integration is not performed on mean boxcar regressors.
%       .fs - sampling frequency of the given input data
%       .HalfDecayTime - time (in s) at which an impulse would be reduced
%           to half its value
% 
%       mirgetdata supercedes get/use_data, but at least one must be
%       specified
% 
% FB 2010.06.18

DEBUG = 1;
dbout = '/home/fbarrett/matlab/logs/fmri_regress_mirfeatures';
check_dir(dbout,0,1);

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
  stim_paths{k} = fileparts(floc{1}{1});
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
  
  % what reg types will we create?
  if isfield(mirinfo.(fname),'reg_types')
    reg_types = mirinfo.(fname).reg_types;
  else
    reg_types = {'raw'};
  end
  nregtype = length(reg_types);
  
  % initialize data matrix & name cell array
  if isfield(mirinfo.(fname),'get_data')
    get_data = mirinfo.(fname).get_data;
    lfnames = {};
    for l=1:length(get_data)
      lfnames = [lfnames fnstr.(get_data{l})];
    end % for l=1:length(get_data)
    lfnames = strcat(fname,lfnames);
  else
    lfnames = {fname};
  end
  nget = length(lfnames);
  funcsig = zeros(npts,nget*nregtype);
  
  % concatenate reg_types and funcnames
  funcnames = cell(1,nget*nregtype);
  for l=1:nregtype
    rni = 1+nget*(l-1):nget*l;
    funcnames(rni) = strcat(lfnames,['_' reg_types{l}]);
  end % for l=1:nregtype
  
  % iterate over stims
  for l=1:nids
    stim_id = stim_ids(l);
    fprintf(1,'feature %s, stim id %d\n',fname,stim_id);

    % search for matching file
    fpath = fullfile(stimulus_root,stim_paths{l},fname);
    fdirs = dir(fullfile(fpath,sprintf('*%d_%s*mat',stim_id,fname)));
    
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
      if compare_structs(job_data.params.mirtoolbox.params,fparams)
        if isfield(job_data,'results') && ~isempty(job_data.results)
          match = 1;
          break
        end
      end
    end % for m=1:length(fdirs
    
    if ~match
      error('could not find %s for stimulus id %d',fname,stim_id);
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
    
    if isfield(mirinfo.(fname),'mirgetdata')
      % initialize the creator function object
      try
        cfh = parse_fh(fdata.data{fcol.creator_func});
        tmp = cfh(1:3);
      end
      
      % load data
      mirdata = load(fdata.data{fcol.mirout_path}{1});
      
      % use mirgetdata to get data
      fsig = mirgetdata(mirdata.a);
    elseif isfield(mirinfo.(fname),'get_data')
      for m=1:nget
        lget = get_data{m};
        ldata = get(mirdata.a,lget);
        if isfield(mirinfo.(fname),'get_proc_fh')
          ldata = get_proc_fh(ldata{1});
        end
        fsig = [fsig ldata];
      end % for m=1:nuse
    end % if isfield(mirinfo.(fname
    
    fsig = repmat(fsig,1,nregtype);
    
    % pre-process all reg types?
    for m=1:nregtype
      rgi = 1+nget*(m-1):nget*m;
      switch reg_types{m}
          case 'raw'
            % leave the signal as-is
          case 'zvar_within_song'
            fsig(:,rgi) = zscore(fsig(:,rgi));
          case 'zvar_across_song'
            % will process later
            error('not yet coded!!!');
          case {'mean','mean_zero'}
            fsig(:,rgi) = repmat(mean(fsig(:,rgi)),size(fsig,1),1);
          otherwise
      end % switch reg_types{m}
    end % for m=1:nregtype

    if DEBUG
      figoutpath = fullfile(dbout,sprintf('%s_%s_%s_stim%d_debug.ps',...
          sess.exp_id,sess.exam_id,fname,stim_id));
      figfh = figure();
      plot(fsig);
      title('raw signal');
      print(figoutpath,'-dpsc','-adobecset');
    end % if DEBUG
    
    % leaky integrator?
    if isfield(minfo,'leaky_integrator') && ~isempty(minfo.leaky_integrator)
      leakfs = minfo.leaky_integrator.fs;
      leakdt = minfo.leaky_integrator.HalfDecayTime;
      li_mask = setxor(1:length(reg_types),strmatch('mean',reg_types));
      fsig(:,li_mask) = IPEMLeakyIntegration(fsig(:,li_mask)',leakfs,leakdt)';
      if DEBUG
        figure(figfh); clf
        plot(fsig);
        title(sprintf('leaky integrated signal: fs %d, dt %d',leakfs,leakdt));
        print(figoutpath,'-append','-dpsc','-adobecset');
      end % if DEBUG
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
    ofset = onset + size(ssig,1) -1;
    funcsig(onset:ofset,:) = ssig;
    
    if DEBUG
      figure(figfh); clf
      plot(ssig);
      title(sprintf('spline interpolated signal: stime %d, TR/dt %d',stime,TR/dt));
      print(figoutpath,'-append','-dpsc','-adobecset');
    end % if DEBUG
  end % for l=1:nids

  % zero-center mean_zero regressors
  zcols = [];
  for l=1:length(funcnames)
    if ~isempty(findstr('mean_zero',funcnames{l})), zcols = [zcols l]; end
  end % for l=1:length(funcnames
  if ~isempty(zcols)
    zrows = funcsig(:,zcols) ~= 0;
    funcsig(zrows,zcols) = zscore(funcsig(zrows,zcols));
    if DEBUG
      figoutpath = fullfile(dbout,sprintf('%s_%s_%s_model%d_debug.ps',...
          sess.exp_id,sess.exam_id,fname,minfo.model_id));
      figure(figfh); clf
      plot(funcsig);
      title('zscored');
      print(figoutpath,'-append','-dpsc','-adobecset');
    end % if DEBUG
  end % if ~isempty(zcols)
  
  % convolve with HRF, add to vals
  convreg = [];
  for l=1:size(funcsig,2)
    convreg(:,l) = conv(funcsig(:,l),spm_hrf(TR/dt));
    vals = [vals convreg([0:actual_nvol-1]*dt+1,l)];
  end % for l=1:size(funcsig,2
  
  if DEBUG
    figoutpath = fullfile(dbout,sprintf('%s_%s_%s_model%d_debug.ps',...
        sess.exp_id,sess.exam_id,fname,minfo.model_id));
    figure(figfh); clf
    plot(vals);
    title('convolved');
    print(figoutpath,'-append','-dpsc','-adobecset');
  end % if DEBUG
  
  names = [names funcnames];
end % for k=1:nmirfunc

% % scale to ~ -1:1
% vals = vals - mean(vals(:));
% vals = vals./max(max(vals(:)),abs(min(vals(:))));
% vals = vals - repmat(mean(vals),size(vals,1),1);
