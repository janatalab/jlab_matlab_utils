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
%           get_data arguments that were used to extract data from a given
%           function output. If mirgetdata is not used here, specifying
%           get_data will direct this script to the data that you want to
%           extract.
%       .(feature function).use_data - cell array of strings identifying
%           those variables in get_data that you want to use. If this is
%           not specified, then use_data defaults to the values in get_data
% 
%       mirgetdata supercedes get/use_data, but at least one must be
%       specified
% 
% FB 2010.06.18

%% init vars
names = {};
vals = [];

m = pinfo.mysql;

TR = pinfo.scanner.TR;
dt = pinfo.scanner.dt;
npts = pinfo.scanner.actual_nvol*dt/TR;

vals = zeros(npts,0);

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
  floc = mysql_extract_data('conn_id', conn_id, ...
      'table','stimulus', ...
      'extract_flds','location', ...
      'stimulus_id',stim_id);
  [stim_paths{k},stim_names{k}] = fileparts(floc{1}{1});
end

% get data from each function for all stims
fprintf(1,'Generating regressors'\n);
for k=1:nmirfunc
  fname = mirfuncs{l};
  fparams = mirinfo.(fname).params;
  
  if isfield(mirinfo.(fname),'get_data')
    get_data = mirinfo.(fname).get_data;
    if ~isfield(mirinfo.(fname),'use_data')
      use_data = get_data;
    else
      use_data = mirinfo.(fname).use_data;
    end
    nuse = length(use_data);
    funcsig = zeros(npts,nuse);
    funcnames = strcat(fname,use_data);
  else
    funcsig = zeros(npts,1);
    funcnames = {fname};
  end
  
  for l=1:nids
    stim_id = stim_ids(l);
    fprintf(1,'feature %s, stim id %d\n',fname,stim_id);

    % search for matching file
    fpath = fullfile(stim_paths{l},fname);
    fdirs = dir(fullfile(fpath,sprintf('*%s*mat',fname)));
    
    match = 0;
    for m=1:length(fdirs)
      if fdirs(m).isdir, continue, end
      clear job_data;
      load(fullfile(fpath,fdirs(m).file));
      if ~exist('job_data','var'), continue, end
      if compare_structs(job_data.params,fparams)
        match = 1;
        break
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
      mirdata = load(fdata.data{fcol.mirout_path});
      
      % use mirgetdata to get data
      fsig = mirgetdata(mirdata);
    elseif isfield(mirinfo.(fname),'get_data')
      for m=1:nuse
        idx = strmatch(use_data{m},get_data);
        fsig = [fsig fdata.data{fcol.get_data}{1}{idx}{1}];
      end % for m=1:nuse
    end % if isfield(mirinfo.(fname

    % spline interpolation!
    stime = mirinfo.(fname).stime;
    svals = 0:stime:(size(fsig,1)-stime);
    ovals = 0:TR/dt:(size(fsig,1)/stime - TR/dt);
    ssig  = [];
    for m=1:size(fsig,2)
      ssig = [ssig spline(svals,fsig(:,m),ovals)];
    end % for m=1:size(fsig
    
    % cut off at music_dur_max?
    if isfield(minfo,'music_dur_max')
      % make sure that curr_val doesn't exceed music_dur_max
      dmax = minfo.music_dur_max*dt/TR;
      if size(ssig,1) > dmax
        ssig(dmax+1:end,:) = [];
      end
    end

    ofset = min(dmax,size(ssig,1));
    funcsig(onsets(l):onsets(l)+ofset,:) = ssig;
  end % for l=1:nids
  vals = [vals funcsig];
  names = [names funcnames];
end % for k=1:nmirfunc

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

% Convolve the regressors with an hrf
for k=1:size(vals,2)
  convreg = conv(vals(:,k),spm_hrf(TR/dt));
  vals(:,k) = convreg(1:size(vals,1));
end

% scale to ~ -1:1
vals = vals - mean(vals(:));
vals = vals./max(max(vals(:)),abs(min(vals(:))));
vals = vals - repmat(mean(vals),size(vals,1),1);
