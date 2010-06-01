function outData = ensemble_summary_subject_stats(inData)
%
% Function that reports summary subject statistics of a given experiment
% such as gender count and ages.
%
% INPUT
%   a cell array of 2 structs: session_info and subject_info, as returned 
%   by ensemble_load_expinfo. This function will be called from within
%   ensemble_load_expinfo, so there shouldn't normally be a need to call
%   this function explicitly. However, in the event that we have data 
%   from an old analysis that didn't include the output from this function, 
%   it can be called easily by passing this information.
%  
% OUTPUT
%   An ensemble data struct with the following variables:
%      nsubs -                    The number of subjects that were included in the analysis
%      mean_age -                 The mean age of included subjects
%      std_age  -                 The standard deviation of ages included in analysis
%      age_range -                The minumum and maximum ages of included subs
%      num_female -               The number of subjects that reported female gender
%      num_male   -               The number of subjects that reported male gender
%      num_gender_no_report -     The number of subjects that didn't report gender
%      sessions_per_sub -         The number of sessions per subject. If the number of
%                                 sessions for all subjects is not equal, the minumum 
%                                 and maximum will be reported.
%      mean_session_dur_minutes - The mean session duration in minutes. If
%                                 there are multiple sessions per subject, 
%                                 the mean for each session order
%                                 will be reported (e.g. if 4 sessions/subject, 
%                                 then 4 values are reported).
%
% Nov 3, 2009 - Stefan Tomic, First Version
% May 8, 2010 - S.T. fixed handling of anon_<hash> subIDs, which don't exist in
%               subject table. Loops by subject ID retrieved in sessInfo
% May 22, 2010  PJ - added searching on type=session_info, and
%               type=subject_info when search for these variables in name field
%               fails. Fixed handling of age=0.
  

fParams.name = 'session_info';
an_idx = ensemble_find_analysis_struct(inData,fParams);

% if searching on name field failed, try searching on type field
if isempty(an_idx)
  fParams = struct('type','session_info');
  an_idx = ensemble_find_analysis_struct(inData,fParams);
end
if isempty(an_idx)
  fprintf('Failed to find session_info\n')
  return
end

sessInfo = inData{an_idx};
sessInfoCols = set_var_col_const(sessInfo.vars);

fParams.name = 'subject_info';
an_idx = ensemble_find_analysis_struct(inData,fParams);
% if searching on name field failed, try searching on type field
if isempty(an_idx)
  fParams = struct('type','subject_info');
  an_idx = ensemble_find_analysis_struct(inData,fParams);
end
if isempty(an_idx)
  fprintf('Failed to find subject_info\n')
  return
end
subInfo = inData{an_idx};
subInfoCols = set_var_col_const(subInfo.vars);

subids = unique(sessInfo.data{sessInfoCols.subject_id});
nsubs = length(subids);
nFemales = 0;
nMales = 0;
nGenderNoReport = 0;

for isub = 1:nsubs
  
  thisSubID  = subids{isub};
  
  [subInSubTable,subInfoIdx] = ismember(thisSubID,subInfo.data{subInfoCols.subject_id});
  
  if(subInSubTable)
    subDOB = subInfo.data{subInfoCols.dob}{subInfoIdx};
  
    if(~isnan(subDOB))
      dobDatenum = datenum(subDOB,'yyyy-mm-dd');
    else
      dobDatenum = NaN;
    end
      
    subGender = subInfo.data{subInfoCols.gender}{subInfoIdx};
    
    switch(subGender)
     case 'F'
      nFemales = nFemales +1;
     case 'M'
      nMales = nMales+1;
     otherwise
      nGenderNoReport = nGenderNoReport + 1;
    end
    
  else
    subDOB = NaN;
    dobDatenum = NaN;
    nGenderNoReport = nGenderNoReport + 1;
  end
    
  sessionIdxs = strmatch(thisSubID,sessInfo.data{sessInfoCols.subject_id},'exact');
  nSess(isub) = length(sessionIdxs);
  
  sessDatenums = sessInfo.data{sessInfoCols.date_time}(sessionIdxs);
  sessEndDatenums = sessInfo.data{sessInfoCols.end_datetime}(sessionIdxs);
  
  
  %use the earliest session that this subject participated in to determine age
  useSessDatenum = min(sessDatenums);
  
  %sort sessions for this sub by start_time to find session order
  [sortedSessDatenums,sortedIdxs] = sort(sessDatenums);
  sortedSessEndDatenums = sessEndDatenums(sortedIdxs);  
  
  nSessThisSub = length(sortedIdxs);
  for iSess = 1:nSessThisSub
    
    serialSessDuration = sortedSessEndDatenums(iSess) - sortedSessDatenums(iSess);
    sessDurations(isub,iSess) = serialSessDuration * 24 * 60;
    
  end
  
  serialAge = useSessDatenum - dobDatenum;
  subAges(isub) = floor(serialAge/365);
  
  % If the age is somehow set to zero, e.g. if subject accidentally entered
  % current day as birthday, enter NaN
  if subAges(isub) == 0
    subAges(isub) = NaN;
  end
  
end

%replace zeros in sessDurations with NaNs. These are sessions that were not
%recorded for a subject (e.g. all subjects completed three sessions, except for
%one subject that only completed two sessions).
sessDurations(sessDurations == 0) = NaN;

%find mean duration per session
nMaxSess = size(sessDurations,2);
for iSess = 1:nMaxSess
  sessMeanDur(iSess) = nanmean(sessDurations(:,iSess));
end

meanAge = nanmean(subAges);
stdAge = nanstd(subAges);
ageRange = [nanmin(subAges) nanmax(subAges)];

if(all(diff(nSess) == 0))
  nSessPerSub = nSess(1);
else
  nSessPerSub = [min(nSess) max(nSess)];
end


outData = ensemble_init_data_struct;
outData.vars = {'nsubs','mean_age','std_age','age_range','num_female','num_male','num_gender_no_report' ...
	   'sessions_per_sub' 'mean_session_dur_minutes'};
outDataCols = set_var_col_const(outData.vars);
outData.data{outDataCols.nsubs} = nsubs;
outData.data{outDataCols.mean_age} = meanAge;
outData.data{outDataCols.std_age} = stdAge;
outData.data{outDataCols.age_range} = ageRange;
outData.data{outDataCols.num_female} = nFemales;
outData.data{outDataCols.num_male} = nMales;
outData.data{outDataCols.num_gender_no_report} = nGenderNoReport;
outData.data{outDataCols.sessions_per_sub} = nSessPerSub;
outData.data{outDataCols.mean_session_dur_minutes} = sessMeanDur;
