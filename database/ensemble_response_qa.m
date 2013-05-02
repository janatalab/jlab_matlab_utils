function out_st = ensemble_response_qa(data_st,params)
% Performs quality-control on response data 
% out_st = ensemble_response_qa(data_st,params);
%
% Calculates the standard deviation of response values across the
% questions on the forms present in the response data that are entered in
% data_st and optionally specified as filtering criteria in the params
% struct.

% 02May2013 Petr Janata

out_st = [];

% Make sure we have a form_id variable
if ~strcmp('form_id',data_st.vars)
  fprintf('%s: No form_id data found. Exiting ...\n', mfilename);
  return
end
cols = set_var_col_const(data_st.vars);

% Check to see whether a minimum standard deviation criterion was
% specified. Otherwise, use zero.
if isfield(params, 'sdCrit')
  sdCrit = params.sdCrit;
else
  sdCrit = 0;
end

% Filter the input data
if isfield(params,'filt')
  data_st = ensemble_filter(data_st, params.filt);
end

% Identify the list of unique forms and make a mask matrix
[formMaskMtx, formids] = make_mask_mtx(data_st.data{cols.form_id}); % Generate a form mask matrix
nforms = length(formids);

% Get a list of subjects
subids = unique(data_st.data{cols.subject_id});
nsubs = length(subids);

% Initialize accumulator variables
out_st.vars = {'subject_ids','form_ids','std_dev','suspect'};
ocols = set_var_col_const(out_st.vars);
out_st.data{ocols.subject_ids} = subids;
out_st.data{ocols.form_ids} = formids;
out_st.data{ocols.std_dev} = nan(nsubs, nforms);


% Loop over subjects
for isub = 1:nsubs
  currSub = subids{isub};
  
  % Create a subject mask
  submask = ismember(data_st.data{cols.subject_id}, currSub);
  
  % Loop over forms
  for iform = 1:nforms
    compositeMask = submask & formMaskMtx(:,iform);
    out_st.data{ocols.std_dev}(isub,iform) = nanstd(enum2data(data_st.data{cols.response_enum}(compositeMask)));
  end % for iform=
end % for isub

% Now check for potentially bad subjects
data = out_st.data{ocols.std_dev};
nanMask = any(isnan(data),2);
sdMask = any(data <= sdCrit, 2);
problemMask = nanMask | sdMask;

out_st.data{ocols.suspect} = subids(problemMask);
nprob = sum(problemMask);
if nprob
  fprintf('%s: Identified %d potentially problematic subjects\n', mfilename, sum(problemMask));

  % Handle subjects with NaNs
  probSubs = subids(nanMask);
  for iprob = 1:sum(nanMask)
    fprintf('Subject %s had no data for form(s): %s\n', probSubs{iprob}, ...
      sprintf('%d, ',formids(isnan(data(strcmp(probSubs{iprob},subids),:)))));
  end
  
  % Handle subjects with no variability
   probSubs = subids(sdMask);
  for iprob = 1:sum(sdMask)
    fprintf('Subject %s had no variability on form(s): %s\n', probSubs{iprob}, ...
      sprintf('%d, ',formids(data(strcmp(probSubs{iprob},subids),:) <= sdCrit)));
  end
 
end

return