function [names,vals] = fmri_regress_tonreg(pinfo,minfo,sess)

% generates tonality tracking regressors
% 
%   [names,vals] = fmri_regress_tonreg(pinfo,minfo,sess)
% 
% 
% 
% FB 2010.01.27

names = {};
vals = [];

regid = pinfo.regid;
m = pinfo.mysql;
gendefs = minfo.ipem;
gendefs.ensemble.conn_id = mysql_make_conn(m);
gendefs.scanner = pinfo.scanner;

% Get a list of the IDs to process
sfilt.include.all = minfo.response_filter;
sdata = ensemble_filter(pinfo,sfilt);
pc = set_var_col_const(pinfo.vars);

proc_ids = cellfun(@str2num,sdata.data{pc.EVENT_CODE});
nids = length(proc_ids);
onsets = sdata.data{pc.RUN_REL_TIME}/1000;

if ~nids
  fprintf(1,'no sound stimuli in this run, SKIPPING\n');
return
end

tonregval = cell(nids,1);

% are we segregating regressors?
if isfield(minfo,'tonreg') && ~isempty(minfo.tonreg)
  % segregation params
  tregp = minfo.tonreg;
  segnames = tregp.seg_names;
  segvals = tregp.seg_vals;
  nseg = length(segvals);
  ccm = parse_fh(minfo.cond_cue_map);
  cue = ccm(tregp.seg_cue);
  resp_params = extract_resp_params_v2(cue,pinfo,minfo,sess);
else
  tregp = '';
end

TR = pinfo.scanner.TR;
dt = pinfo.scanner.dt;  % actually number of steps rather than delta time
npts_detail = pinfo.scanner.actual_nvol*pinfo.scanner.dt;

for j = 1:nids
  fprintf('Generating tonality regressors for stim %d/%d\n',j,nids);

  if proc_ids(j) > length(tonregval)
    tonregval{j} = [];
  end

  [tonregval{j},tonregnames] = fmri_generate_tonreg(proc_ids(j),gendefs);
  curr_val = tonregval{j};
  
  if isfield(minfo,'music_dur_max')
    % make sure that curr_val doesn't exceed music_dur_max
    dmax = minfo.music_dur_max*dt/TR;
    if size(curr_val,1) > dmax
      curr_val(dmax+1:end,:) = [];
    end
  end

  if j == 1
    curr_val_cols = size(curr_val,2);
    if isstruct(tregp)
      num_tonreg = curr_val_cols*nseg;
    else
      num_tonreg = curr_val_cols;
    end
    tonregmtx = zeros(npts_detail,num_tonreg);
  end

  start_idx = round(onsets(j)*dt/TR);
  stop_idx = start_idx+size(curr_val,1)-1;

  % Determine which columns we are putting these tonality regressors
  % into. This applies if we are breaking up regressors by
  % autobiographical memory.
  curr_cols = [];
  if isstruct(tregp)
    rval = resp_params(j);
    seg = 0;
    for k=1:nseg
      if any(ismember(segvals{k},rval))
        seg = k;
        break
      end
    end
    if seg
      curr_cols = curr_val_cols*(seg-1)+1:curr_val_cols*seg;
    end
  else
    curr_cols = 1:num_tonreg;
  end
  if ~isempty(curr_cols)
    tonregmtx(start_idx:stop_idx,curr_cols) = curr_val;
  end
end % for j = 1:nids

% Convolve the regressors with an hrf
convreg = [];
for ireg = 1:size(tonregmtx,2)
  convreg(:,ireg) = ...
      conv(tonregmtx(:,ireg),spm_hrf(TR/dt));
  if isstruct(tregp)
    curr_name = sprintf('%s_%s',...
        tonregnames{mod(ireg-1,length(tonregnames))+1},...
        segnames{fix((ireg-1)/length(tonregnames))+1});
  else
    curr_name = tonregnames{ireg};
  end
  names{ireg} = curr_name;
  vals(:,ireg) = convreg([0:(pinfo.scanner.actual_nvol-1)]*dt+1,ireg);
end

% % scale to ~ -1:1
% vals = vals - mean(vals(:));
% vals = vals./max(max(vals(:)),abs(min(vals(:))));
% vals = vals - repmat(mean(vals),size(vals,1),1);

% zscore song timepoints
vals(vals ~= 0) = zscore(vals(vals ~= 0));
