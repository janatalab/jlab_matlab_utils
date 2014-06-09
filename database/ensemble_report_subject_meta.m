function out_st = ensemble_report_subject_meta(data_st,params)
% Summarizes subject information, e.g. age, gender
%
% Requires subject_info structure such as that returned by
% mysql_get_subinfo or mysql_get_expinfo, as well as a session_info
% structure. Both of the required structures are obtained most easilty by
% mysql_get_expinfo().
%

% 31Dec2012 Petr Janata
% 09Jun2014 PJ - added option to filter, and added returning of data struct

out_st = ensemble_init_data_struct;
out_st.vars = {'subject_id','age','gender'};
ocols = set_var_col_const(out_st.vars);
out_st.data = cell(1,length(out_st.vars));

% Find the subject info analysis 
idx = ensemble_find_analysis_struct(data_st, struct('type','subject_info'));
if isempty(idx)
  error('Could not find subject info structure');
end

subst = data_st{idx};
subcols = set_var_col_const(subst.vars);

% Perform any filtering on the subject_info structure
if isfield(params,'filt')
  subst = ensemble_filter(subst, params.filt);
end

% Find the session info structure
idx = ensemble_find_analysis_struct(data_st, struct('type','session_info'));
if isempty(idx)
  error('Could not find session info structure')
end
sessst = data_st{idx};
sesscols = set_var_col_const(sessst.vars);

% Perform any filtering on the session_info structure
if isfield(params,'filt')
  sessst = ensemble_filter(sessst, params.filt);
end
fid = 1;

%
% Determine the age at each session. 
%

% Match the subject ID associated with each session with the DOB entry in
% the subject info structure
[~, srcidx] = ismember(sessst.data{sesscols.subject_id}, subst.data{subcols.subject_id});

expdate = sessst.data{sesscols.date_time};
dob = subst.data{subcols.dob}(srcidx);
female_mask = strcmp(subst.data{subcols.gender}(srcidx),'F');
agevect_s = etime(datevec(expdate),datevec(dob));
agevect_y = agevect_s/(60*60*24*365);
fprintf(fid,'N=%d (%d females)\n', length(agevect_y), sum(female_mask));
fprintf(fid,'Age range: %d - %d\n', fix(min(agevect_y)), fix(max(agevect_y)));
fprintf(fid,'Age (mean +/- std): %2.1f (%2.1f)\n', mean(agevect_y), std(agevect_y));

% Populate the output structuresubst.data{subcols.gender}(srcidx)
out_st.data{ocols.subject_id} = subst.data{subcols.subject_id}(srcidx);
out_st.data{ocols.gender} = subst.data{subcols.gender}(srcidx);
out_st.data{ocols.age} = agevect_y;


return