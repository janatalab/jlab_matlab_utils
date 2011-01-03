function durs = fmri_stim_duration(pinfo,minfo,sids)

% returns durations for given stimuli
% 
%   durs = fmri_stim_duration(sids)
% 
% REQUIRES
%   pinfo.mysql.conn_id - open database connection id
%   minfo.music_dur_max (optional)
%   sids - vector of stimulus ids
% 
% FB 2010.02.25

durs = [];

% get durations from the database
local_conn_id = 0;
try conn_id = pinfo.mysql.conn_id;
catch
    conn_id = 0;
    local_conn_id = 1;
    mysql_make_conn([],[],conn_id);
end

dstr = sprintf('SELECT duration FROM stimuli WHERE stimulus_id IN (%s)',...
  regexprep(num2str(sids'),'\s\s',','));
times = mysql(conn_id,dstr);
if local_conn_id
    mysql(conn_id,'close');
end

if length(times) ~= length(sids)
    error('wrong number of durations returned from the database')
end
durs = str2num(datestr(times,'HH'))*360+...
  str2num(datestr(times,'MM'))*60+str2num(datestr(times,'SS'));

if isfield(minfo,'music_dur_max')
    % set upper limit to durations
    durs(durs > minfo.music_dur_max) = minfo.music_dur_max;
end
