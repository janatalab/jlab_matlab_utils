function analysis_list = ensemble_jobman(analysis_list,params)
% Top-level job manager for the Ensemble Matlab analysis toolkit.
%
% analysis_list = ensemble_jobman(analysis_list,params)
%
% analysis_list is a cell array of analysis structures that
% conform to the ensemble analysis structure spec (see wiki for details).
% possible params are:
%
% params.ensemble.conn_id: integer that specifies a mysql connection to use
%                          the jobman will use this id to open a connection
% params.run_analyses:     a vector of indices that specify which analysis
%                          jobs to run. Analysis jobs whose indices are not
%                          listed will not run, but their results (from a
%                          previous run) will be available.
%
% params.resultNameIsAnalysisName: if True, results.name will be set to the
%                          name of the analysis
%
%
%
% This function serves as a controller for the Ensemble Matlab
% analysis toolkit. It processes an analysis_list and calls the
% appropriate analysis functions as specified by the function
% handles, and with the appropriate parameters specified in the
% params struct.Results from the analysis are stored in the
% 'results' field of the analysis list structure.
%


% 16 Feb 2007 - First version, built on code by PJ, Stefan Tomic
% 01 May 2007 - PJ Added better host, database, conn_id support
% 11/27/07 - improved conn_id support
% 8/25/08  - ST, if conn_id is empty, conn_local is now set to
%            false. Reverted last version, since 'dbstop if error' didn't take
%            you to the desired breakpoint
% 06/15/10 - PJ cleaning mysql_make_conn support
% 04/11/11 - PJ added alternative loading of results struct from file
% 10Nov2012 - PJ added and fixed insertion of analysis name into the
%             results structure
% 02Sep2014 - Fixed handling of run_analyses vector passed in as column
%             rather than row vector
% 01Sep2016 - Added option to force the name of the result structure to
%             take on the analysis name


%if conn_id is set in params, then open a mysql connection
%this mysql connection may be used by analysis functions that
%specify the same conn_id in their params
conn_local = false;

try
  conn_id = params.ensemble.conn_id;
catch
  conn_id = [];
  conn_local = false;
end

try ignoreEmptyConnID = params.ensemble.ignoreEmptyConnID;
catch ignoreEmptyConnID = 0;
end

% If connection ID is not empty, check to see if it is open
if ~isempty(conn_id) && ~ignoreEmptyConnID
  if mysql(conn_id,'status')
    if conn_id == 0
      conn_local = true;
    end
    mysql_make_conn(params.ensemble);
  end
end

if(isfield(params,'run_analyses'))
  idxs = params.run_analyses;
else
  idxs = 1:length(analysis_list);
end

for aidx = 1:length(idxs)
  ia = idxs(aidx);
  analysis_name = analysis_list{ia}.name;
  
  fprintf('\nPerforming analysis %d/%d: %s\n', ia, length(analysis_list), analysis_name);
  indata = {};
  
  %assign the input data structure for the function based on the
  %required fieldname
  if isfield(analysis_list{ia},'requires') && ~isempty(analysis_list{ia}.requires)
    
    for requiredIdx = 1:length(analysis_list{ia}.requires)
      requiredData = analysis_list{ia}.requires{requiredIdx};
      anSearchCrit.name = requiredData.name;
      anListIdx = ensemble_find_analysis_struct(analysis_list,anSearchCrit);
      if(isempty(anListIdx))
        error(sprintf('Analysis %s specifies a non-existent required analysis %s',...
          analysis_list{ia}.name,anSearchCrit.name));
      end
      
      numIdx = length(anListIdx);
      if numIdx > 1
        error('Found %d analysis structs (%s) that matched the name of the required structure: %s', ...
          numIdx, sprintf('%d ', anListIdx), anSearchCrit.name)
      end
      
      if(isfield(analysis_list{anListIdx},'results') )
        %search for the required data type
        
        dataList = analysis_list{anListIdx}.results;
        
        % Handle the possibility that the results is simply a flag and that the
        % results are stored to disk, e.g. if the analysis was run using
        % ensemble_jobman_parallel.
        if isnumeric(dataList) && dataList == 1 && ...
            isfield(analysis_list{anListIdx},'path') && ~isempty(analysis_list{anListIdx}.path)
          tmp = load(analysis_list{anListIdx}.path);
          dataList = tmp.job_data.results;
        end
        
        %see if 'vars' was specified, in which case we need to
        %search for the variable to return
        if(isfield(requiredData,'vars') && ~isempty(requiredData.vars))
          
          
          %find the index of the result we need
          dataListIdx = strmatch(requiredData.vars,dataList.vars);
          
          %if result data struct was found then assign to indata
          %otherwise, throw an error indicating the data is not there.
          if(~isempty(dataListIdx))
            indata{requiredIdx} = dataList.data{dataListIdx};
            
          else
            error([sprintf('Analysis %s specifies a non-existent required',analysis_name)...
              sprintf(' datatype %s for analysis %s',requiredData.vars,analysis_list{anListIdx}.name)]);
          end
          
        else
          %if the 'vars' field of the required data list was not
          %specified, then return the whole results data structure
          indata{requiredIdx} = analysis_list{anListIdx}.results;
        end
        
        
      else
        
        %results field wasn't populated so required data wasn't available
        error(sprintf('Required data for analysis ''%s'' was not available.',analysis_name));
        
      end
      
    end % for requiredIdx
    
  end %if isfield(analysis_list{ia},'requires')
  
  %if there was only one required data struct, then don't pass this
  %as a cell
  if(length(indata) == 1)
    indata = indata{1};
  end
  
  % Call the function that is registered with this analysis
  fh = analysis_list{ia}.fun;
  
  if(isfield(analysis_list{ia},'params'))
    analysisParams = analysis_list{ia}.params;
  else
    analysisParams = [];
  end
  
  result = fh(indata,analysisParams);
  
  % Set the name to the analysis name if isfield(analysisParams,'resultNameIsAnalysisName')
  if isfield(analysisParams,'resultNameIsAnalysisName') && analysisParams.resultNameIsAnalysisName
    result.name = analysis_name;
  end
  
  analysis_list{ia}.results = result;
  
  % Check to see if the name of the result is empty, and if so, populate it
  % with the name of this analysis
  numResults = length(result);
  for ires = 1:numResults
    if isstruct(result(ires))
      if ~isfield(result(ires),'name') || isempty(result(ires).name) || isfield(params,'resultNameIsAnalysisName') && params.resultNameIsAnalysisName
        analysis_list{ia}.results(ires).name = analysis_name;
      else
        continue
      end
    elseif ~isfield(analysis_list{ia}.results{ires},'name') || isempty(analysis_list{ia}.results{ires}.name) || isfield('params','resultNameIsAnalysisName') && params.resultNameIsAnalysisName
      analysis_list{ia}.results{ires}.name = analysis_name;
    end
  end
  %    check_conn_exit(params)
  %    err = lasterror;
  %    fprintf('Failed on analysis %d ...\n%s\nReturning completed analyses ...\n', ia, err.message);
  %    return
  
end % for ia - analysis list

if conn_local
  mysql(conn_id,'close');
end
end
