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

out_st = ensemble_init_data_struct;
out_st.type = 'response_time_data';
out_st.vars = [{'data_transform','rank'} params.rt.calc];

% Get the idx of the response_data
ridx = ensemble_find_analysis_struct(data_st,struct('type','response_data'));
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

% Precalculate the masks we are using for this looping variable
if strcmp(item_type,'none')
  item_mask_mtx = ones(length(rst.data{1}),1);
  itemids = 1;
else
  % Precalculate the masks for the item variable we are using
  [item_mask_mtx, itemids] = make_mask_mtx(rst.data{rcols.(item_str)});
end
nitems = length(itemids);

% Convert the response text data to doubles
allRTs = cellfun(@str2num,rst.data{rcols.response_text});

% Convert to seconds
allRTs = allRTs/1000;

% Now loop over items
for iitem = 1:nitems
  item_mask = item_mask_mtx(:,iitem);
  
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
  [ranked(ixfm).mean, ranks] = sort(out_st.data{ocols.mean}(xfm_mask),'descend');
  
  % Write the rank info to the by-item data matrix
  out_st.data{ocols.rank}(xfm_mask,1) = ranks(:); 
  
  for ivar = 1:length(rankVars)
    cv = rankVars{ivar};
    
    % Rank the current variable
    tmp = out_st.data{ocols.(cv)}(xfm_mask);
    ranked(ixfm).(cv) = tmp(ranks);
  end

  % Calculate the standard error of the mean
  ranked(ixfm).sem = ranked(ixfm).std./sqrt(ranked(ixfm).numel-1);
  
  % Generate a rank ordered bar graph
  subplot(nxfm,1,ixfm)
  h = bar(ranked(ixfm).mean);
  add_errorbars(h,ranked(ixfm).sem)
  
  % Add xtick labels
  set(gca,'xtick',1:nitems,'xticklabel','')
  
  for iitem = 1:nitems
    if iscell(ranked(ixfm).(item_str))
      itemStr = sprintf('%s',ranked(ixfm).(item_str){iitem});
    else
      itemStr = sprintf('%d',ranked(ixfm).(item_str)(iitem));
    end
    text(iitem,-1,itemStr,'rotation',-90)
  end
  
  xlabel(strrep(item_str,'_','\_'))
  
end

end % function

