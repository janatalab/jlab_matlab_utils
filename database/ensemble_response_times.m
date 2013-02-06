function [out_st] = ensemble_response_times(data_st,params)
% Analyzes response time values stored in the response_text field of a
% response table.
%
% [out_st] = ensemble_response_times(data_st,params);
%
% The questions for which response times are to be analyzed should be
% specified in a filtering directive stored in the params structure, e.g.
% filt.include.all.question_id = [question_ids];
%
% params.rt.rt_fld = name of field containing response time data (in msec).
% Default is response_text
%
% params.rt.multiRespCheckFld = name of field to use for checking
% uniqueness of response to each stimulus, i.e. one response per stimulus
% per subject. Default='session_id';
%
% params.rt.transform = {'none','log10'}; - specifies whether a transform
% should be applied to the RT data. Default='none'.
%
% params.rt.proc = {}; - list of functions to process. The default behavior is to
% treat the name as a function handle. Default processing includes
% calculation of the mean and standard deviation.
%
% params.rt.by_item = {'stimulus_id','subject_id'};
%
% NOTE: A cleaner implementation of this would require that data are
% already pre-filtered. What is required for that to happen though is a
% modification of ensemble_filter() selectively remove rows based on
% combinations of variable values.


% 18Nov2012 Petr Janata
% 05Feb2013 PJ - added flexible specification of field containing response
%                times

out_st = ensemble_init_data_struct;
out_st.type = 'response_time_data';
out_st.vars = [{'data_transform','rank'} params.rt.calc];

% Get the idx of the response_data
ridx = ensemble_find_analysis_struct(data_st,struct('type','response_data'));
if isempty(ridx)
  ridx = ensemble_find_analysis_struct(data_st,struct('type','resp_x_stim'));
end
if isempty(ridx)
  error('Could not find an appropriate data structure with response data');
end
rst = data_st{ridx};
rcols = set_var_col_const(rst.vars);

% Apply an standard filtering parameters
if isfield(params,'filt')
  fprintf('Applying filtering criteria\n')
  rst = ensemble_filter(rst, params.filt);
end

% See if there is a datastruct containing filtering parameters and apply
% them
if isfield(params.rt,'exclude')
  rm_vars = {'name','type','report','meta'};
  idx = ensemble_find_analysis_struct(data_st,struct('type','filtering_data'));
  
  fparams = rmfield(ensemble_tree2datastruct(data_st{idx}), rm_vars);
  nfilt = size(fparams.data{1},1);
  
  % Have to filter row by row
  for ifilt = 1:nfilt
    tmp.vars = fparams.vars;
    for ivar = 1:length(fparams.vars)
      tmp.data{ivar} = fparams.data{ivar}(ifilt);
    end
    f.exclude.all = rmfield(ensemble_datastruct2tree(tmp), rm_vars);
    rst = ensemble_filter(rst,f);
  end
end

% Determine whether we are looping by a particular item (stimulus_id,
% subject_id)
if isfield(params,'item_var')
  item_str = params.item_var;
else
  item_str = 'stimulus_id';
end

switch item_str
  case 'session_id'
    item_type = 'session';
  case 'subject_id'
    item_type = 'subject';
  case 'stimulus_id'
    item_type = 'stimulus';
  case 'none'
    item_type = 'none';
  otherwise
    fprintf('%s: Unknown item variable: %s\n', mfilename, item_str);
end

% Set the output structure name
out_st.name = sprintf('by_%s', item_str);

% Add item types to the list of output variables
out_st.vars = [item_str out_st.vars];
ocols = set_var_col_const(out_st.vars);

% If we are looping over stimulus_id and we are normalizing the response
% time data by some value, make sure we have a structure of those values
if strcmp(item_type,'stimulus') && isfield(params.rt,'transform') && any(ismember(params.rt.transform,'norm'))
  % Get the stimulus info
  sidx = ensemble_find_analysis_struct(data_st,struct('type','stimulus_metadata'));
  if isempty(sidx)
    error('stimulus_metadata structure required for analysis, but none found\n')
  end
  stim_st = data_st{sidx};
  stimcols = set_var_col_const(stim_st.vars);
end

% Precalculate the masks we are using for this looping variable
if strcmp(item_type,'none')
  item_mask_mtx = ones(length(rst.data{1}),1);
  itemids = 1;
else
  % Precalculate the masks for the item variable we are using
  [item_mask_mtx, itemids] = make_mask_mtx(rst.data{rcols.(item_str)});
end
nitems = length(itemids);

% Determine which field the response times are in
if isfield(params.rt, 'rt_fld')
  rt_fld = params.rt.rt_fld;
else
  rt_fld = 'response_text';
end

% Convert the response text data to doubles
if iscell(rst.data{rcols.(rt_fld)})
  allRTs = cellfun(@str2num,rst.data{rcols.(rt_fld)});
else
  allRTs = rst.data{rcols.(rt_fld)};
end

% Convert to seconds
allRTs = allRTs/1000;

if isfield(params.rt, 'multiRespCheckFld')
  multiRespCheckFld = params.rt.multiRespCheckFld;
else
  multiRespCheckFld = 'session_id';
end

% Now loop over items
for iitem = 1:nitems
  item_mask = item_mask_mtx(:,iitem);
  
  % Check to see if we have multiple responses per subject per item
  uniqueSess = unique(rst.data{rcols.(multiRespCheckFld)}(item_mask));
  if iscell(uniqueSess)
    [~,cnt] = cellhist(uniqueSess);
  else
    [cnt] = hist(rst.data{rcols.(multiRespCheckFld)}(item_mask),uniqueSess);
  end
  
  if any(cnt > 1)
    multmask = ismember(rst.data{rcols.(multiRespCheckFld)}, uniqueSess(cnt > 1));
    uniqueMult = unique(rst.data{rcols.(multiRespCheckFld)}(multmask));
    nmult = length(uniqueMult);
    for imult = 1:nmult
      currMult = uniqueMult(imult);
      currmask = ismember(rst.data{rcols.(multiRespCheckFld)}, currMult) & item_mask;

      fprintf('Found %d responses for stimulus %d, session %d\n', ...
        sum(currmask), itemids(iitem), currMult);

      % Take the first of the multiple responses
      [~,minidx] = min(rst.data{rcols.date_time}(currmask));
      curridxs = find(currmask);
      goodmask = false(size(currmask));
      goodmask(curridxs(minidx)) = true;
      item_mask = xor(item_mask,xor(currmask,goodmask));
    end
  end
  
  % Extract the response times
  resptimes = allRTs(item_mask);
   
  % Determine whether we are going to transform the data
  if isfield(params.rt,'transform')
    xfms = params.rt.transform;
    nxfm = length(xfms);
  else
    nxfm = 1;
    xfms = {'none'};
  end
  
  % Loop over possible transformations of the data
  for ixfm = 1:nxfm
    xfmType = xfms{ixfm};
    
    row = (iitem-1)*nxfm + ixfm; % Determine which row we are putting data into
    
    % Output information about the item ID
    if iscell(itemids)
      out_st.data{ocols.(item_str)}{row,1} = itemids(iitem);
    else
      out_st.data{ocols.(item_str)}(row,1) = itemids(iitem);
    end
    
    % Output information about the transform type
    out_st.data{ocols.data_transform}{row,1} = xfmType;
    
    % Transform the data
    switch xfmType
      case 'none'
        data = resptimes;
      case 'norm' % normalize by item-level data
        % Normalize by a particular field
        if strcmp(params.rt.norm,'duration') && isnan(stim_st.data{stimcols.duration}(iitem))
          srcmask = stim_st.data{stimcols.stimulus_id} == itemids(iitem);
          stimname = fullfile(params.paths.stimulus_root, stim_st.data{stimcols.location}{srcmask});
          system_str = sprintf('mp3info -p %%S %s', stimname);
          [~,r] = system(system_str);
          stim_st.data{stimcols.duration}(srcmask) = min(params.maxListenTime,str2double(r));
        end

        normVal = stim_st.data{stimcols.(params.rt.norm)}(stim_st.data{stimcols.stimulus_id} == itemids(iitem));
        data = resptimes/normVal;
      otherwise
        fh = str2func(xfmType);
        data = fh(resptimes);
    end
    
    % Calculate the desired metrics on these response times
    ncalc = length(params.rt.calc);
    for icalc = 1:ncalc
      currCalc = params.rt.calc{icalc};
      
      switch currCalc
        
        otherwise
          fh = str2func(currCalc);
          out_st.data{ocols.(currCalc)}(row,1) = fh(data);
          
      end % switch currCalc
    end % for icalc
  end % for ixfm
end % for iitem

% Perform additional analyses across items by transform type
rankVars = {'numel','std',item_str};
for ixfm = 1:nxfm
  xfmType = xfms{ixfm};
  xfm_mask = ismember(out_st.data{ocols.data_transform}, xfmType);
  
  % Get the rank order of the mean data
  [ranked(ixfm).mean, srcidxs] = sort(out_st.data{ocols.mean}(xfm_mask),'descend');
  
  % Write the rank info to the by-item data matrix
  [sorted_srcidx, ranks] = sort(srcidxs);
  out_st.data{ocols.rank}(xfm_mask,1) = ranks; 
  
  for ivar = 1:length(rankVars)
    cv = rankVars{ivar};
    
    % Rank the current variable
    tmp = out_st.data{ocols.(cv)}(xfm_mask);
    ranked(ixfm).(cv) = tmp(srcidxs);
  end

  % Calculate the standard error of the mean
  ranked(ixfm).sem = ranked(ixfm).std./sqrt(ranked(ixfm).numel-1);
  
  % Generate a rank ordered bar graph
  subplot(nxfm,1,ixfm)
  h = bar(ranked(ixfm).mean,'facecolor',ones(1,3)*0.8);
  add_errorbars(h,ranked(ixfm).sem,'k')
  
  % Add xtick labels
  set(gca,'xtick',1:nitems,'xticklabel','')
  
  for iitem = 1:nitems
    if iscell(ranked(ixfm).(item_str))
      itemStr = sprintf('%s',ranked(ixfm).(item_str){iitem});
    else
      itemStr = sprintf('%d',ranked(ixfm).(item_str)(iitem));
    end
    text(iitem,max(get(gca,'ylim'))*0.05,itemStr,'rotation',90)
  end
  
  xlabel(strrep(item_str,'_','\_'), 'fontsize', 14, 'fontweight', 'bold')
  
  switch lower(xfmType)
    case 'none'
      lblstr = 'Response time (s)';
    otherwise
      lblstr = sprintf('%s(s)', xfmType);
  end
  ylabel(lblstr, 'fontsize', 14, 'fontweight','bold')
  
  % Add sample size info
  if ~any(diff(ranked(ixfm).numel))
    nobsstr = sprintf('N = %d',ranked(ixfm).numel(1));
  else
    nobsstr = sprintf('N = %d - %d', min(ranked(ixfm).numel), max(ranked(ixfm).numel));
  end
  text(0.95,0.95,nobsstr,'units','normalized',...
    'horizontalalign','right', ...
    'fontsize', 16, ...
    'fontweight','bold')
end

end % function

