function ensemble_copy_stimuli2dir(varargin)
% Copies stimuli specified either in data_st or in the experiment specified
% in params to a desired output directory.
%
% ensemble_copy_stimuli2dir(varargin)
%
% Takes a series of tag/value pairs. The following tags are recognized:
%
% REQUIRES either
% {'experiment_title', 'experiment_name', 'experiment'} - the value in the
%      experiment_title field of the experiment table
% OR
% {'stimulus_id'} - an array of stimulus IDs to export
% OR
% {'response_table'} - name of response table to pull stimulus IDs from
% {'conn_id'} - REQUIRED - an active connection ID to the database
%
% OPTIONAL
% {'filt','filter'} - a structure compatible with ensemble_filter() 
% {'outpath'} - path where stimuli should be copied to. Default is
%      './stimuli'
%
% This script is useful for sharing stimuli used in an experiment with
% other people.

% 09Jan2014 Petr Janata

% Initialize required arguments
expname = '';
resptbl = '';
stimids = [];
conn_id = [];
outroot = './stimuli';
stimroot = '';
filt_st = [];

% Process input arguments
for iarg = 1:2:nargin
  currArg = varargin{iarg};
  switch currArg
    case {'experiment','experiment_name','experiment_title'}
      expname = varargin{iarg+1};
      
    case {'conn_id'}
      conn_id = varargin{iarg+1};
      
    case {'stimulus_id','stimids'}
      stimids = varargin{iarg+1};
      
    case {'response_table'}
      resptbl = varargin{iarg+1};
      
    case {'filt','filter'}
      filt_st = varargin{iarg+1};
      
    case {'stimroot'}
      stimroot = varargin{iarg+1};
      
    case {'outpath','outdir','outroot'}
      outroot = varargin{iarg+1};
      
    otherwise
      fprintf('Input argument (%s) not recognized\n', currArg);
  end % switch currArg
end % for iarg

if isempty(stimids) && isempty(resptbl) && isempty(expname)
  error('Must specify one of the following: stimulus_id, response_table, experiment_name')
end

if isempty(conn_id) || mysql_check_conn(conn_id)
  error('Must specify conn_id for MySQL database connection. The conn_id must be active')
end

if isempty(stimroot)
  error('Must specify root location of stimuli (stimroot)')
end

% We need a list of stimulus IDs
if isempty(stimids)
  % We need a response table name
  if isempty(resptbl)
    fprintf('Fetching response_table name from experiment: %s\n', expname);
    mysql_str = sprintf('SELECT response_table FROM experiment WHERE experiment_title="%s";', expname);
    resptbl = mysql(conn_id, mysql_str);
    resptbl = resptbl{1};
  end % if isempty(resptbl)
  
  % Get the stimuli from the response table
  fprintf('Fetching distinct stimulus IDs from response_table: %s\n', resptbl);
  mysql_str = sprintf('SELECT DISTINCT stimulus_id FROM %s;', resptbl);
  stimids = mysql(conn_id, mysql_str);
  
end % if isempty(stimids)

% Eliminate any NaNs from the stimids and make sure they are unique
stimids = unique(stimids(~isnan(stimids)));

% Get the stimulus info
fprintf('Extracting stimulus metadata for %d stimuli\n', length(stimids));
stimMeta = mysql_extract_metadata('table','stimulus', ...
  'stimulus_id',stimids, ...
  'conn_id', conn_id);

% Convert the stimMeta tree to an Ensemble data structure
stim_st = ensemble_tree2datastruct(stimMeta,struct('ignore_reserved_names',true));

% Perform any desired filtering
if ~isempty(filt_st)
  stim_st = ensemble_filter(stim_st, filt_st);
end

scols = set_var_col_const(stim_st.vars);

% Make sure the output directory exists
check_dir(outroot);

% Copy the stimuli
srcs = strcat(fullfile(stimroot,filesep), stim_st.data{scols.location});
nstim = size(srcs,1);
for istim = 1:nstim
  srcLoc = srcs{istim}; % get the current stimulus path
  
  [~,destStub,ext] = fileparts(srcLoc); % get the current stimulus name
  
  destLoc = fullfile(outroot,[destStub ext]); % figure out where it's going
  
  % Copy it
  unix_str = sprintf('cp %s %s', srcLoc, destLoc);
  status = unix(unix_str);
  if status
    error('Problem executing: %s', unix_str)
  end
  
end % for istim
  
return