function PL = set_col_const(var_list)
% Figures out which data columns in a Presentation logfile correspond to
% different variables 
%
% var_list is a cell string containing field names
% PL is a structure with constants (column ids) assigned for the following
% fields:
%
% SUB_ID -- subject ID
% STIM_ID -- stimulus ID
% DATE_TIME -- datestamp
% QUEST_ID -- question ID
% QUEST_TXT -- question text
% SUBQUEST_TXT -- subquestion text (heading)
% SUBQUEST_ID -- subquestion ID
% QUEST_ITER -- question iteration
% RESP_ENUM -- response value
% RESP_TXT -- response text
% STIM_TYPE -- stimulus type
% PAIR_IDX -- pair index
% 

% 04/11/05  Petr Janata - adapted from the set_col_const.m script I wrote for
%                         database stuff
% 05/03/11 - FB - added STIM_TYPE, PAIR_IDX (new output vars)

% The mappings variable has structure field names in the 1st column and
% variable names in the 2nd column.  Abandon the mappings approach for the time being.
mappings = { ...
      'SUB_ID', 'Subject'; ...
      'TRIAL_NUM', 'Trial'; ...
      'EVENT_TYPE', 'Event Type'; ...
      };


PL = struct('SUB_ID',[], ...
    'TRIAL_NUM',[], ...
    'EVENT_TYPE',[], ...
    'EVENT_CODE',[], ...
    'EVENT_ABSTIME',[], ...
    'EVENT_TTIME',[], ...
    'TTIME_UNC',[], ...
    'DUR_UNC',[], ...
    'REQ_TIME',[], ...
    'REQ_DUR',[], ...
    'RUN_REL_TIME',[], ...
    'RESP_CODE',[], ...
    'RESP_TIME',[], ...
    'RUN',[], ...
    'STIM_ID',[], ...
    'STIM_TYPE',[], ...
    'PAIR_IDX',[] ...
    );

nvars = length(var_list);
for ivar = 1:nvars
  switch var_list{ivar}
    case {'Subject','SUB_ID'}
      PL.SUB_ID = ivar;
      
    case {'Trial','TRIAL_NUM'}
      PL.TRIAL_NUM = ivar;
  
    case {'Event Type','EVENT_TYPE'}
      PL.EVENT_TYPE = ivar;
      
    case {'Code','EVENT_CODE'}
      PL.EVENT_CODE = ivar;
    
    case {'Time','EVENT_ABSTIME'}
      PL.EVENT_ABSTIME = ivar;  % absolute time in experiment - in msec (10^-3)
      
    case {'TTime','EVENT_TTIME'}
      PL.EVENT_TTIME = ivar;  % time in trial - 10^-4

    case {'Stim Type','STIM_TYPE'}
      PL.STIM_TYPE = ivar;
      
    case {'Pair Index','PAIR_IDX'}
      PL.PAIR_IDX = ivar;
          
    case 'Uncertainty'
      switch var_list{ivar-1}
	case {'TTime','TTIME_UNC'}
	  PL.TTIME_UNC = ivar;
	case {'Duration','DUR_UNC'}
	  PL.DUR_UNC = ivar;
	otherwise
	  fprintf('set_col_const: Uncertain about uncertainty column\n');
      end
    
    case {'ReqTime','REQ_TIME'}
      PL.REQ_TIME = ivar; % 10^-4
    
    case {'ReqDur','REQ_DUR'}
      PL.REQ_DUR = ivar; % 10^-4
      
      case {'ResponseCode','RESP_CODE'}
          PL.RESP_CODE = ivar;
          
      case {'ResponseTime','RESP_TIME'}
          PL.RESP_TIME = ivar;
          
      case {'Run Relative Time','RUN_REL_TIME'}
          PL.RUN_REL_TIME = ivar;
          
      case {'Run','RUN'}
          PL.RUN = ivar;

      case {'stimulus id','stimulus_id','stim_id','STIM_ID'}
          PL.STIM_ID = ivar;
          
  end % switch
end % for ivar

return