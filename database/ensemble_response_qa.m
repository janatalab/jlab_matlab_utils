function out_st = ensemble_response_qa(data_st,params)
% Performs quality-control on response data 
% out_st = ensemble_response_qa(data_st,params);
%
% Calculates the standard deviation of response values across the
% questions on the forms present in the response data that are entered in
% data_st and optionally specified as filtering criteria in the params
% struct.

% 02May2013 Petr Janata
% 05Jul2013 PJ - added formatting for convenient output string
% 21Oct2015 PJ - enhanced analyses and reporting - optional check for
%                minimum number of iterations of specified forms

out_st = [];

% Make sure we have a form_id variable
if ~strcmp('form_id',data_st.vars)
  fprintf('%s: No form_id data found. Exiting ...\n', mfilename);
  return
end
cols = set_var_col_const(data_st.vars);

% Set up reporting
if isfield(params,'report')
  fid = ensemble_init_fid(params.report);
else
  fid = 1;
end

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
out_st.vars = {'subject_ids','form_ids','std_dev','suspect','insufficient','numCleanIter'};
ocols = set_var_col_const(out_st.vars);
out_st.data{ocols.subject_ids} = subids;
out_st.data{ocols.form_ids} = formids;
out_st.data{ocols.std_dev} = nan(nsubs, nforms);

if isfield(params,'form')
  minimumIterationForms = fieldnames(params.form);
else
  minimumIterationForms = {};
end
if ~isempty(minimumIterationForms)
  for iform = 1:length(minimumIterationForms)
    minIterFormIDs(iform) = params.form_id_const.(minimumIterationForms{iform});
  end
else
  minIterFormIDs = [];
end
nminIter = length(minIterFormIDs);
out_st.data{ocols.insufficient} = nan(nsubs, nminIter);
out_st.data{ocols.numCleanIter} = nan(nsubs, nminIter);

% Get an overall nanmask
nanmask = isnan(data_st.data{cols.response_enum});

% Loop over subjects
for isub = 1:nsubs
  currSub = subids{isub};
  
  % Create a subject mask
  submask = ismember(data_st.data{cols.subject_id}, currSub);
  
  % Loop over forms
  for iform = 1:nforms
    compositeMask = submask & formMaskMtx(:,iform);
    out_st.data{ocols.std_dev}(isub,iform) = nanstd(enum2data(data_st.data{cols.response_enum}(compositeMask)));
    
    % Check whether this form is on the list of minimum iteration checks
    fidx = find(minIterFormIDs == formids(iform));
    if ~isempty(fidx)
      minCleanIter = params.form.(minimumIterationForms{fidx}).minimumCleanIterations;
      numCleanIter = length(unique(data_st.data{cols.response_order}(compositeMask & ~nanmask)));
      out_st.data{ocols.numCleanIter}(isub,fidx) = numCleanIter;
      
      if numCleanIter < minCleanIter
        out_st.data{ocols.insufficient}(isub,fidx) = 1;
      end
    end
    
  end % for iform=
end % for isub

% Now check for potentially bad subjects
data = out_st.data{ocols.std_dev};
nanMask = any(isnan(data),2);
sdMask = any(data <= sdCrit, 2);
insuffMask = any(out_st.data{ocols.insufficient},2);
problemMask = nanMask | sdMask | insuffMask;

out_st.data{ocols.suspect} = subids(problemMask);
nprob = sum(problemMask);
if nprob
  fprintf(fid,'%s: Identified %d potentially problematic subjects\n', mfilename, sum(problemMask));

  % Handle subjects with NaNs
  probSubs = subids(nanMask);
  fprintf(fid, '\n\n%d subjects missing responses on some forms\n', length(probSubs));
  for iprob = 1:sum(nanMask)
    fprintf(fid, 'Subject %s missing a total of %d responses in encounters with form(s): %s\n', probSubs{iprob}, ...
      sum(isnan(data(strcmp(probSubs{iprob}, subids),:))), ...
      sprintf('%d, ',formids(isnan(data(strcmp(probSubs{iprob},subids),:)))));
  end
  outstr = sprintf('''%s'',', probSubs{:});
  outstr(end) = '';
  fprintf(fid, '\ndefs.suspect.missingResponses = {%s};', outstr);
  
  % Handle subjects with no variability
  probSubs = subids(sdMask);
  fprintf(fid, '\n\n%d subjects with no variability on some forms\n', length(probSubs));
  for iprob = 1:sum(sdMask)
    fprintf(fid, 'Subject %s had no variability on form(s): %s\n', probSubs{iprob}, ...
      sprintf('%d, ',formids(data(strcmp(probSubs{iprob},subids),:) <= sdCrit)));
  end
  outstr = sprintf('''%s'',', probSubs{:});
  outstr(end) = '';
  fprintf(fid, '\ndefs.suspect.responseVariability = {%s};', outstr);
    
  % Handle subjects with an insufficient number of trials
  probSubs = subids(insuffMask);
  fprintf(fid, '\n\n%d subjects with fewer than %d iterations of required forms\n', length(probSubs), minCleanIter);
  for iprob = 1:sum(insuffMask)
    fprintf(fid, 'Subject %s had insufficient (< %d) iterations of form(s): %s\n', probSubs{iprob}, ...
      minCleanIter, ...
      cell2str(minimumIterationForms(out_st.data{ocols.insufficient}(strcmp(subids,probSubs{iprob}),:)),','));
  end
  outstr = sprintf('''%s'',', probSubs{:});
  outstr(end) = '';
  fprintf(fid, '\ndefs.suspect.tooFewTrials = {%s};', outstr); 
  
end

if fid > 1
  fclose(fid);
end


return