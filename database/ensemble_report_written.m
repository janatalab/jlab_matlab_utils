function an_st = ensemble_report_written(data_st,params)
% Displays written responses associated text or varchar fields.
%
% outdata = ensemble_report_written(data_st,params);
%
% Displays written responses associated text or varchar fields.  Currently, the
% reporting is rather primitive, but should be enhanced to support various
% tabular formats.  Table format could be specified by providing a vector of
% compqids to appear in the columns.  A structure array could actually provide
% qid info as well as formatting info, e.g. column width, etc.
%
% Note: Currently, the script will not treat the same question appearing on
% different forms as a different instance of the question. If you don't want
% answers to the same question on different forms combined, you must filter the
% data to only process forms with unique question IDs.  This behavior may
% change in future versions.
%
% params.ensemble.conn_id - mysql connection to use
% params.filt - any filtering criteria that are to be applied
% params.report - structure containing directives regarding the output.
%                 See ensemble_init_fid()
% params.report.print - generate output
% params.report.write2file - write output to a file (otherwise stdout)
% params.report.fname - name of file to write responses to
% params.report.filemode - defaults to 'wt'
% params.report.verbose - defaults to 1

% 04/29/07 Petr Janata - adapted from ensemble_enum_stats
% 07/17/08 PJ - fixed call to ensemble_init_fid to conform to new param handling
% 23Nov2011 PJ - improved handling of conn_id

try
	if isfield(params, 'mysql')
		conn_id = params.mysql.conn_id;
	elseif isfield(params, 'ensemble')
		conn_id = params.ensemble.conn_id;
	end
catch ME
	conn_id = [];
end

an_st = ensemble_init_data_struct;
an_st.type = 'enum_written_by_compqid'; 

% Make sure we have a compqid variable
data_st = ensemble_check_compqid(data_st);
if isempty(data_st)
  return
end

% Set the column constants
incol = set_var_col_const(data_st.vars);

% Apply any specified filtering to the input data
if isfield(params,'filt')
  fprintf('Applying filtering criteria\n')
  data_st = ensemble_filter(data_st, params.filt);
end

%
% Gather metadata on the questions
%

% Get a list of unique composite question IDs
fprintf('Getting list of unique composite question IDs\n')
%qids = fix(unique(data_st.data{incol.compqid}));
qids = unique(fix(data_st.data{incol.compqid}));

fprintf('Extracting question metadata\n');
qinfo = mysql_extract_metadata('table','question', ...
    'question_id',qids, ...
    'conn_id', conn_id);

% Figure out which of the questions in the qinfo structure are varchar or text and remove
% those that are not
fprintf('Removing non-text questions\n');
qinfo_char_mask = ismember({qinfo.type},{'varchar','text'});
if sum(~qinfo_char_mask)
  qinfo(~qinfo_char_mask) = [];
end

% Create masks for all of the response data
char_compqids = [qinfo.compqid];
  
% Filter the data again
filt.include.any.compqid = char_compqids;
data_st = ensemble_filter(data_st,filt);

% Precalculate the subject masks
%[sub_mask_mtx, subids] = make_mask_mtx(data_st.data{incol.subject_id});
%nsub = length(subids);

%
% Check to see if we are writing all responses to a single file and open a file
% if necessary
%
if ~isfield(params,'report')
  params.report = struct();
end
composite_fid = ensemble_init_fid(params.report);

% See if we are only writing responses to the file and not question or subject
% information
try response_only = params.report.response_only; catch response_only = false; end

%
% Loop over all of the unique question/subquestion combinations or compqids
%
nqid = length(qinfo);
for iqid = 1:nqid
  curr_id = qinfo(iqid).compqid;
  an_st.vars{iqid} = sprintf('compqid %s',num2str(curr_id));
  
  % Copy the question info over to the display parameter structure in case we
  % are going to display some of the data
  params.display.qinfo = qinfo(iqid);

  % Filter the data based on the current ID
  clear qfilt
  qfilt.include.all.compqid = curr_id;
  curr_st = ensemble_filter(data_st,qfilt);
  col = set_var_col_const(curr_st.vars);

  if ~response_only
    fprintf(composite_fid,'\n\nCompqid: %.2f\nQUESTION: %s\nSUBQUESTION: %s\n', ...
      curr_id, qinfo(iqid).question_text, qinfo(iqid).heading);
  end
  
  % Loop over the rows in the data and write out each response to the file
  nrows = length(curr_st.data{col.subject_id});
  for irow = 1:nrows
    if ~response_only
      subject_str = sprintf('%s (%d): ', ...
        curr_st.data{col.subject_id}{irow}, ...
        curr_st.data{col.session_id}(irow));
    else
      subject_str = '';
    end
    fprintf(composite_fid,'\n%s%s\n', subject_str, curr_st.data{col.response_text}{irow});
  end
end % for iqid

if composite_fid > 1
  fclose(composite_fid);
end

an_st.meta.params = params;

end % function ensemble_report_written

