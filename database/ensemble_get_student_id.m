function outdata = ensemble_get_student_id(indata,params)

% returns unique student id/subject id combinations found in response data
% 
%   outdata = ensemble_get_student_id(indata,params)
% 
% given an ensemble data struct containing response table data, this
% function will find all unique combinations of subject_id and student_id,
% and will return these in an ensemble data struct.
% 
% REQURIES
%   indata - ensemble data struct containing response table data (see
%       ensemble_load_expinfo)
%   params.student_id_qid (default: 558) - question id containing a
%       participant's student_id
% 
% FB 2009.05.16

% init vars
outdata = ensemble_init_data_struct();
outdata.name = 'student_ids';
outdata.type = 'student_ids';
outdata.vars = {'subject_id','student_id'};
outdata.data = {{} {}};
ocol = set_var_col_const(outdata.vars);

% check structure of indata
incol = set_var_col_const(indata.vars);
if ~isfield(incol,'subject_id'), error('no subid found within indata'); end
if ~isfield(incol,'response_text'), error('response text not found'); end

% get student_ids
if isfield(params,'student_id_qid') && ~isempty(params.student_id_qid)
  siqid =  params.student_id_qid;
else
  siqid =  558;
end
filt.include.all.question_id = siqid
indata = ensemble_filter(indata,filt);

% check for rows
if isempty(indata.data{1})
  error('no student IDs (question %d) found in the given dataset',siqid);
end

% find unique subjects
[sm,us] = make_mask_mtx(indata.data{incol.subject_id});
ns = length(us);
if any(sum(sm) > 1), warning('multiple records for some subjects'); end

% iterate over subjects, find student_id
for is=1:ns
  subid = us{is};
  sids = indata.data{incol.response_text}(sm(:,is));
  if length(sids) > 1
    warning('multiple student ids found for %s, taking first',subid);
    sids = sids(1);
  end
  outdata = ensemble_add_data_struct_row(outdata,'subject_id',subid,...
      'student_id',sids);
end
