function outData = check_completion_v2(varargin)
% provides details on subjects who have started and who have completed a given experiment
%
% adapted from anps_collect_check_completion.m written by Fred Barrett
%
% REQUIRED PARAMS
% Either 
%    params.resptbl_name              - The name of the response table for the experiment
% or 
%    params.ensemble.experiment_title - The title of the experiment
%
% OPTIONAL PARAMS
% params.terminal_form             - The form ID to check against for
%                                    completion of experiment
% params.terminal_question         - The question ID to check against for
%                                    completion
% params.filt                      - Filter params to filter out unwanted subjects
% params.ensemble.conn_id          - Connection ID to use for database 
%                                    if no live connection, a new connection is established
% params.ensemble.host             - hostname of database
% params.ensemble.database         - database name
% params.report_after              - any sessions before this date will be
%                                    ignored in the report (if not specified,
%                                    all sessions will display)
%
% OUTPUTS
%
% Returns a list of subjects who completed the experiment
%
% Author(s)
% 2/21/2009 Fred Barrett - Original Author
% 4/7/2009 Stefan Tomic - adapted script into a more general purpose function for all experiments

outData = [];

if nargin == 1
  params = varargin{1};
elseif nargin == 2
  params = varargin{2};
else
  fprintf('%s: Wrong number of arguments: %d\n', mfilename, nargin);
end

ensemble_globals;
enc_key = ensemble_get_encryption_key;

try
  conn_id = params.ensemble.conn_id;
catch
  conn_id = 7;
  params.ensemble.conn_id = conn_id;
end

try
  host = params.ensemble.host;
catch
  host = '';
end

try
  database = params.ensemble.database;
catch
  database = '';
end

try
  params.filt;
catch
  params.filt = {};
end

if(mysql(conn_id,'status') ~= 0)
  mysql_make_conn(host,database,conn_id);
end

try tbl_name = params.resptbl_name; catch tbl_name = ''; end

if isempty(tbl_name) && isfield(params.ensemble, 'experiment_title')
  mysql_str = sprintf('SELECT response_table, experiment_id FROM experiment WHERE experiment_title="%s"', params.ensemble.experiment_title);
  [tbl_name, exp_id] = mysql(conn_id, mysql_str);
  tbl_name = tbl_name{1};
end

try form_id = params.terminal_form; catch form_id = []; end;

% If form_id is unspecified, get it from the experiment_x_form table
if isempty(form_id)
  mysql_str = sprintf('SELECT form_id, form_order FROM experiment_x_form WHERE experiment_id=%d;', exp_id);
  [ids,order] = mysql(conn_id, mysql_str);
  % Get next to last form, as last form is an end_session form that doesn't
  % leave its mark on the response table
  form_id = ids(order==(max(order)-1));
end

try question_id = params.terminal_question; catch question_id = []; end
if isempty(question_id)
  mysql_str = sprintf('SELECT question_id, form_question_num FROM form_x_question WHERE form_id=%d;', form_id);
  [qids, order] = mysql(conn_id,mysql_str);
  question_id = qids(order==max(order));
end

try PRINT_TO_FILE = params.report.write2file; catch PRINT_TO_FILE=0; end

% % % REPORT AFTER THIS DATE
try
  report_after = datenum(params.report_after);
catch
  %if report_after wasn't specied, report all sessions
  report_after = 0;
end

%% Consult the session table to get list of completed subs, sessions, and time stamps
vars = {'session_id','date_time','end_datetime','subject_id','ticket_id'};
mysql_str = sprintf('SELECT %s FROM session WHERE experiment_id =%d;', cell2str(vars,','), exp_id);
data = cell(1,length(vars));
[data{:}] = mysql(conn_id,mysql_str);
completedData = ensemble_init_data_struct;
completedData.vars = vars;
completedData.data = data;
cdCols = set_var_col_const(vars);

% Get the ticket codes

%% Filter completion data if desired
completedData = ensemble_filter(completedData,params.filt);

% open logfile
if exist('PRINT_TO_FILE','var') && PRINT_TO_FILE
  logfname = fullfile(params.paths.logpath,'completion_info.txt');
  fid = fopen(logfname,'wt');
  fprintf(fid,'Completion information for experiment: %s\n\n\n', params.ensemble.experiment_title);
else
  fid = 1; % stdout
end
params.report.fid = fid;

nsess = length(completedData.data{cdCols.session_id});
fprintf(fid,'%d initiated sessions\n', nsess);

completed_mask = ~isnan(completedData.data{cdCols.end_datetime});
completed_idxs = find(completed_mask);
fprintf(fid,'%d completed sessions\n', sum(completed_mask));

for itype = 1:2
  if itype == 1
    idx_list = completed_idxs;
    type_str = 'COMPLETED';
  else
    idx_list = find(~completed_mask);
    type_str = 'INCOMPLETE';
  end
  
  fprintf(fid,'\n%s\n', type_str);
  fprintf(fid,'Session\tSubject\tTicketID\tTicketCode\n');
  for isess = 1:length(idx_list)
    curr_idx = idx_list(isess);
    
    curr_tick_id = completedData.data{cdCols.ticket_id}(curr_idx);
    mysql_str = sprintf('SELECT ticket_code FROM ticket WHERE ticket_id = %d;', curr_tick_id);
    ticket_code = mysql(conn_id, mysql_str);
    
    fprintf(fid, '%d\t%s\t%d\t%s\n', ...
      completedData.data{cdCols.session_id}(curr_idx), ...
      completedData.data{cdCols.subject_id}{curr_idx}, ...
      curr_tick_id, ticket_code{1});
  end
end


return
%% EVERYTHING BELOW HERE NEED TO BE FIXED


%% Print out list of those session and subject IDs that completed
curr_filt.exclude.any.end_datetime = NaN;
ensemble_display_table(ensemble_filter(completedData,curr_filt),params.report);

[completedSubs,completedSess,completedDateTimes] = ...
    mysql(conn_id,sprintf(['select subject_id,session_id,date_time from' ...
		    ' %s where form_id = %d and question_id = %d group by session_id'],tbl_name,form_id,question_id));

%place in ensemble data struct so that we can filter out the subs
%that we aren't interested in.
completedData = ensemble_init_data_struct;
completedData.vars = {'subject_id','session_id','date_time'};
completedDataCols = set_var_col_const(completedData.vars);
completedData.data{completedDataCols.subject_id} = completedSubs;
completedData.data{completedDataCols.session_id} = completedSess;
completedData.data{completedDataCols.date_time} = completedDateTimes;
completedData = ensemble_filter(completedData,params.filt);

subids = completedData.data{completedDataCols.subject_id};
sessids = completedData.data{completedDataCols.session_id};
timestamps = completedData.data{completedDataCols.date_time};
nsess = length(subids);


fprintf(fid,'participants who have completed\n');
fprintf(fid,'-------------------------------\n');
fprintf(fid,'Name\t\t\tSubID\t\tSession\tTimestamp\t\tNumStims\tDuration\n');
for isess = 1:nsess
    if strmatch('tmp_',subids{isess}), continue, end
    idx = strmatch(subids{isess},subs.data{sCol.subject_id},'exact');
    if timestamps(isess) > report_after
        if isempty(subs.data{sCol.name_last}(idx)) || ...
                isempty(subs.data{sCol.name_last}(idx)) || ...
                isempty(subs.data{sCol.name_first}(idx))
          fprintf(fid,'Cannot find any subjects\n');
        end
        subid = subids{isess};
        first = cell2str(subs.data{sCol.name_first}(idx));
        last  = cell2str(subs.data{sCol.name_last}(idx));
        nstim = subs.data{sCol.nstim}(idx);
        dur   = subs.data{sCol.dur}(idx)*24;
        fprintf(fid,'%s %s\t\t%s\t%d\t%s\t%d\t\t%.2f\n',...
		first, last,subids{isess},sessids(isess),...
		datestr(timestamps(isess)),nstim,dur);
    end
end

fprintf(fid,'\n\nTOTAL COMPLETED:%d\n\n',isess);

ndidxs  = ~ismember(subs.data{sCol.subject_id},subids);
notdone = subs.data{sCol.subject_id}(ndidxs);

% Tally # completed stim and total time spent for all ppts
expinfo = ensemble_load_expinfo({},params);
expCols = set_var_col_const(expinfo.vars);
subs = expinfo.data{expCols.subject_info};
resp = expinfo.data{expCols.response_data};
sCol = set_var_col_const(subs.vars);
rCol = set_var_col_const(resp.vars);
newcols = {'nstim','dur'};
for inc = 1:length(newcols);
    nidx = length(fieldnames(sCol))+1;
    subs.vars{nidx} = newcols{inc};
    subs.data{nidx} = [];
    sCol.(newcols{inc}) = nidx;
end
nsub = length(subs.data{sCol.subject_id});

stimFilt.exclude.any.stimulus_id = NaN;

fprintf(fid,'participants who have started\n');
fprintf(fid,'-----------------------------\n');
fprintf(fid,'Name\t\t\tSubID\t\tNumStims\tDuration\n');
for isub=1:nsub
    subid = subs.data{sCol.subject_id}(isub);
    subFilt.include.all.subject_id = subid;
    subResp = ensemble_filter(resp,subFilt);
    if (~isempty(subResp.data{1}))
        start = min(subResp.data{rCol.date_time});
        stop  = max(subResp.data{rCol.date_time});
	stimResp = ensemble_filter(subResp,stimFilt);
        sstim = unique(stimResp.data{rCol.stimulus_id});
    else
        start = 0;
        stop  = 0;
        sstim = [0];
    end
    subs.data{sCol.dur}(isub)   = stop-start;
    subs.data{sCol.nstim}(isub) = length(unique(sstim));
    if isempty(subs.data{sCol.name_last}(isub)) || ...
            isempty(subs.data{sCol.name_last}(isub)) || ...
            isempty(subs.data{sCol.name_first}(isub)) || ...
            isempty(subid) || isempty(subs.data{sCol.nstim}(isub)) || ...
            isempty(subs.data{sCol.dur}(isub))
      error('empty item\n');
    end
    fprintf(fid,'%s\t\t%s\t%d\t\t%1.2f\n',sprintf('%s %s',...
        cell2str(subs.data{sCol.name_first}(isub)),...
        cell2str(subs.data{sCol.name_last}(isub))),...
        cell2str(subid),...
        subs.data{sCol.nstim}(isub),(subs.data{sCol.dur}(isub)*24));
end

fprintf(fid,'\n\n');


fprintf(1,'DONE!\n\n');

return

% Get subject information
subid_str = sprintf('''%s'',', subids{:});
subid_str(end) = [];

mysql_str = sprintf(['SELECT s.subject_id, aes_decrypt(`s`.`name_last`,''' enc_key '''), ' ...
		    'aes_decrypt(`s`.`name_first`,''' enc_key '''), ' ...
    'count(DISTINCT(r.stimulus_id)) FROM subject s, %s r '...
    'WHERE s.subject_id = r.subject_id AND s.subject_id IN (%s) '...
    'GROUP BY r.subject_id'], tbl_name, subid_str);

[unique_subids, last, first, stim] = mysql(conn_id, mysql_str);

% Print out the information
nsess = length(unique_subids);

fprintf(fid,'Name\tSubID\tSession\tTimestamp\tStim Completed\tisess\n');
for isess = 1:nsess
    idx = strmatch(subids{isess}, unique_subids,'exact');
    if timestamps(isess) > report_after
        fprintf(fid,'%s\t%s\t%d\t%s\t%d\t%d\n',sprintf('%s %s',first{idx},last{idx}),...
          subids{isess},sessids(isess),datestr(timestamps(isess)),stim(isess),isess);
    end
end

fprintf(fid,'\n\nTOTAL COMPLETED:%d\n\n',isess);

% % % %  those that may have not completed the final form

excl_str = cell2str(params.filt.exclude.any.subject_id,''',''');

mysql_str = ['select aes_decrypt(`s`.`name_last`,''' enc_key '''), ' ...
	     'aes_decrypt(`s`.`name_first`,''' enc_key '''),' ... 
	     'r.subject_id,r.session_id,' ...
	     'count(distinct(r.stimulus_id)) from ' tbl_name ' r, subject s '...
	     'where s.subject_id = r.subject_id and r.subject_id not like (''' ...
	     excl_str ''') group by r.subject_id'];

[chklast,chkfirst,chksubid,chksess,chkstim] = mysql(conn_id,mysql_str);

nchk = length(chklast);

n31 = 0;
working.first = {};
working.last  = {};
working.subid = {};
working.sess  = [];
working.stim  = [];

fprintf(fid,'Name\tSubid\tSessIDs\tStim Completed\n');
for ichk = 1:nchk
    if ((chkstim(ichk) > 24) && ...
            ~length(strmatch(cell2str(chksubid(ichk)),subids)))
        fprintf(fid,'%30s %s\t%s\t%s\t%s\n',cell2str(chkfirst(ichk)),...
            cell2str(chklast(ichk)),cell2str(chksubid(ichk)),...
            num2str(chksess(ichk)),num2str(chkstim(ichk)));
            n31 = n31+1;
    elseif(chkstim(ichk) < 30)
        working.first = [working.first chkfirst(ichk)];
        working.last  = [working.last  chklast(ichk)];
        working.subid = [working.subid chksubid(ichk)];
        working.sess  = [working.sess  chksess(ichk)];
        working.stim  = [working.stim  chkstim(ichk)];
    end
end

fprintf(fid,'\n\nTOTAL STARTED:%d\n\nTOTAL not caught:%d\n\n',nchk,n31);

for iwk = 1:length(working.first)
   fprintf(fid,'%30s %s\t%s\t%d\t%d\n',cell2str(working.first(iwk)),...
       cell2str(working.last(iwk)),cell2str(working.subid(iwk)),...
       working.sess(iwk),working.stim(iwk));
end

fprintf(fid,'\n\nTOTAL IN PROGRESS:%d\n\n',length(working.first));

% fprintf(fid,'\n\nTOTAL COMPLETED:%d\n\n',isess)
% 
% fprintf(fid,'\n\nTOTAL STARTED:%d\n\nTOTAL not caught:%d\n\n',nchk,n31)
