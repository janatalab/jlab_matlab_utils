function outData = check_completion_v2(params)
% provides details on subjects who have started and who have completed a given experiment
%
% adapted from anps_collect_check_completion.m written by Fred Barrett
%
% REQUIRED PARAMS
% params.resptbl_name              - The name of the response table for the experiment
% params.terminal_form             - The form ID to check against for
%                                    completion of experiment
% params.terminal_question         - The question ID to check against for completion
% params.filt                      - Filter params to filter out unwanted subjects
% params.ensemble.experiment_title - The title of the experiment

% OPTIONAL PARAMS
% params.ensemble.conn_id          - Connection ID to use for database 
%                                    if no live connection, a new connection is established
% params.ensemble.host             - hostname of database
% params.ensemble.database         - database name
%
%
% OUTPUTS
%
% Returns a list of subjects who completed the experiment
%
% Author(s)
% 2/21/2009 Fred Barrett - Original Author
% 4/7/2009 Stefan Tomic - adapted script into a more general purpose function for all experiments

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

if(~mysql(conn_id,'status'))
  mysql_make_conn(host,database,conn_id);
end

tbl_name = params.resptbl_name;

form_id = params.terminal_form;
question_id = params.terminal_question;

PRINT_TO_FILE=0;

% % % REPORT AFTER THIS DATE
% report_after = datenum('01-Dec-2008 12:00:00');
report_after = datenum(params.report_after);

% open logfile
if exist('PRINT_TO_FILE','var') && PRINT_TO_FILE
  logfname = fullfile(params.paths.logpath,...
      sprintf('anps_cc_%s.txt',datestr(now,30)));
  fid = fopen(logfname,'wt');
else
  fid = 1; % stdout
end

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

%get list of completed subs, sessions, and time stamps
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
          error('empty item\n');
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

fprintf(1,'DONE!\n\n');

return

% Get subject information
subid_str = sprintf('''%s'',', subids{:});
subid_str(end) = [];

% SELECT s.subject_id, s.name_last, s.name_first, count(DISTINCT(r.stimulus_id))
% FROM subject s, response_gpref r
% WHERE s.subject_id = r.subject_id
% AND s.subject_id IN (%)
% GROUP BY r.subject_id

mysql_str = sprintf(['SELECT s.subject_id, s.name_last, s.name_first, '...
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

% select s.name_last, s.name_first, r.subject_id,count(distinct(r.stimulus_id)) 
% from response_gpref r, subject s
% where s.subject_id = r.subject_id
% and s.subject_id not like '01ttf%'
% group by r.subject_id

excl_str = cell2str(params.filt.exclude.any.subject_id,''',''');

mysql_str = ['select s.name_last,s.name_first,r.subject_id,r.session_id,' ...
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
