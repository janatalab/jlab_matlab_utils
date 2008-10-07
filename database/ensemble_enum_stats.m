function an_st = ensemble_enum_stats(data_st,params)
% Calculates statistics on responses to enum questions.
% 
% outdata = ensemble_enum_stats(data_st,params);
%
% Calculates various descriptive and quantitative statistics on responses to
% questions that are enums.
%
% Note: Currently, the script will not treat the same question appearing on
% different forms as a different instance of the question. If you don't want
% answers to the same question on different forms combined, you must filter the
% data to only process forms with unique question IDs.  This behavior may
% change in future versions.
%
% 02/02/07 Petr Janata - adapted from ensemble_enum_hist

an_st = ensemble_init_data_struct;
an_st.type = 'enum_stats_by_compqid'; 

% Make sure that a stats structure has been specified as part of the params structure
if ~isfield(params,'stats')
  fprintf('ensemble_enum_stats: No stats requests specified\n');
  return
end

% Extract info regarding the database we should be talking to
try database = params.ensemble.database; catch database = 'ensemble_main'; end

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
qids = fix(unique(data_st.data{incol.compqid}));

qinfo = mysql_extract_metadata('database', database, 'table','question','question_id',qids);

% Figure out which of the questions in the qinfo structure are enums and remove
% those that are not
qinfo_enum_mask = ismember({qinfo.type},'enum');
if sum(~qinfo_enum_mask)
  fprintf('ensemble_enum_stats: Removing %d non-enum qids\n', sum(~qinfo_enum_mask));
  qinfo(~qinfo_enum_mask) = [];
end

% Figure out which of the questions are bitmasks that allow for selection of
% multiple values (checkbox as opposed to radiogroup) and remove these from the
% qinfo array
qinfo_bitmask_mask = ismember({qinfo.html_field_type},'checkbox');
if sum(qinfo_bitmask_mask)
  fprintf('ensemble_enum_stats: Removing %d bitmask qids from list of qids\n', sum(qinfo_bitmask_mask));
  qinfo(qinfo_bitmask_mask) = [];
end

% Create masks for all of the response data
enum_compqids = [qinfo.compqid];
  
% Filter the data again
filt.include.any.compqid = enum_compqids;
data_st = ensemble_filter(data_st,filt);

% Copy all of the enum data that are not a bitmask to the data_vect and convert
% to category indices
data_vect = data_st.data{incol.response_enum};
data_vect = enum2data(data_vect);

% Precalculate the subject masks
[sub_mask_mtx, subids] = make_mask_mtx(data_st.data{incol.subject_id});
nsub = length(subids);

%
% Set stuff up for writing to a file, if that's what we're going to do.
%
%fid = ensemble_init_fid(params.display.tables);

%
% Loop over all of the unique question/subquestion combinations or compqids
%
nqid = length(qinfo);
for iqid = 1:nqid
  an_st.vars{iqid} = sprintf('compqid %s',num2str(qinfo(iqid).compqid));
  
  % Copy the question info over to the display parameter structure in case we
  % are going to display some of the data
  params.display.qinfo = qinfo(iqid);

  % Get the enum categories
  enum_values = qinfo(iqid).enum_values;
  ncat = length(enum_values);
  
  % Make a mask for the data corresponding to this question
  qid_mask = ismember(data_st.data{incol.compqid},qinfo(iqid).compqid);
    
  an_st_l1 = ensemble_init_data_struct;
  an_st_l1.type = 'enum_basic_stats';
  an_st_l1.vars = {'by_subject','across_subjects'};
  an_st_l1.meta.question = qinfo(iqid);

  an_st_l1_cols = set_var_col_const(an_st_l1.vars);
  nlevel1 = length(an_st_l1.vars);
  
  tmp_st = {};
  tmp_idx = [];
  for il1 = 1:nlevel1
    id_str = an_st_l1.vars{il1};
    tmp_idx.(id_str) = il1;
    
    tmp_st{il1} = ensemble_init_data_struct;
    tmp_st{il1}.type = sprintf('enum_stats_%s', id_str);

    % The variables in this analysis are all of the fields within the stats
    % field of the params structure
    stats_list = fieldnames(params.stats.(id_str));
    nstats = length(stats_list);
    
    switch id_str
      case 'by_subject'
	aux_vars = {'subject_id','nresp'};
      case 'across_subjects'
	aux_vars = {'nsub'};
    end
    tmp_vars = [aux_vars stats_list'];
    tmp_st_cols = set_var_col_const(tmp_vars);
    tmp_st{il1}.vars = tmp_vars;
  
    % Initialize output variables
    for ia = 1:length(tmp_vars)
      switch tmp_vars{ia}
	case {'subject_id'}
	  tmp_st{il1}.data{ia} = subids;
	otherwise
	  switch id_str
	    case 'by_subject'
	      tmp_st{il1}.data{ia} = zeros(nsub,1);
	    otherwise
	      tmp_st{il1}.data{ia} = [];
	  end
      end
    end % for ia=1:length(tmp_vars)

    %
    % Execute some type-specific code. One can imagine additional case
    % statements for 'by_trial' or 'by_attribute'
    %
    
    switch id_str
      %
      % Deal with the set of by_subject analyses
      %
      case 'by_subject'
	for isub = 1:nsub
	  sub_mask = sub_mask_mtx(:,isub);

	  % Tally the number of responses the subject made to this question
	  nresp = sum(sub_mask&qid_mask);
	  tmp_st{il1}.data{tmp_st_cols.nresp}(isub) = nresp;
	  
	  if ~nresp
	    no_resps = 1;
	  else
	    no_resps = 0;
	  end
	  
	  % Now loop over all of the analyses we want to perform
	  for istat = 1:nstats
	    stat_str = stats_list{istat};
	    
	    % If we need to enter a Nan, do that here
	    if no_resps
	      tmp_st{il1}.data{tmp_st_cols.(stat_str)}(isub) = NaN;
	      continue
	    end
	    
	    switch stat_str
	      case {'mean','std','min','max'}
		fh = str2func(stat_str);
		tmp_st{il1}.data{tmp_st_cols.(stat_str)}(isub) = fh(data_vect(sub_mask&qid_mask));
	    end % switch stat_str
	  end % for istat=
	end % for isub
	
      case 'across_subjects'
	src_st = tmp_st{tmp_idx.by_subject};
	src_cols = set_var_col_const(src_st.vars);

	% 02/02/07 PJ Currently hard-coded to use subject-level means as input into this
	% level of the analysis. Ultimately, this should really become another
	% level of abstraction which supports different types of source data.
	src_data = src_st.data{src_cols.mean};
	      
	% Remove any data with NaNs
	src_data(any(isnan(src_data),2),:) = [];
	
	tmp_st{il1}.data{tmp_st_cols.nsub} = size(src_data,1);
	for istat = 1:nstats
	  stat_str = stats_list{istat};
	  switch stat_str
	    case {'mean','std','min','max'}
	      
	      % Evaluate the basic function
	      fh = str2func(stat_str);
	      tmp_st{il1}.data{tmp_st_cols.(stat_str)} = fh(src_data);
	      
	      % See if there is additional processing to be done
	      if isstruct(params.stats.(id_str).(stat_str))
		proc_list = fieldnames(params.stats.(id_str).(stat_str));
		for iproc = 1:length(proc_list)
		  switch proc_list{iproc}
		    case 'ttest'
		      try mu = params.stats.(id_str).(stat_str).ttest.mu; ...
			  catch mu = 'midpoint'; end
		      if isstr(mu) && strcmp(mu,'midpoint')
			mu = (ncat+1)/2;
		      end
		      
		      tmp_st2 = ensemble_init_data_struct;
		      tmp_st2.type = proc_list{iproc};
		      tmp_st2.vars = {'H','p','ci','stats'};
		      [tmp_st2.data{1:nargout(proc_list{iproc})}] = ttest(src_data, mu);
		    otherwise
		      continue
		  end
		  tmp_st{il1}.vars{end+1} = sprintf('%s_%s',stat_str,proc_list{iproc});
		  tmp_st_cols = set_var_col_const(tmp_st{il1}.vars);
		  tmp_st{il1}.data{end+1} = tmp_st2;
		end % for iproc
	      end % if isstruct(params.stats.(id_str).(stat_str)
	  end % switch stat_str
	end % for istat=
		
      otherwise
	
    end % switch id_str (by_subject, across_subjects)
    
    % Register a reporting function and execute it if desired
    tmp_st{il1}.report.fun = str2func(sprintf('report_stats_%s', id_str));
    
    try do_report = params.report.print_tables; catch do_report = 1; end
    if do_report
      fprintf('Doing report for %s\n', func2str(tmp_st{il1}.report.fun))
      params.report.question = qinfo(iqid);  % cludge
      tmp_st{il1}.report.fun(tmp_st{il1},params.report);
    end
  end % for il1 = 1:nlevel1
  an_st_l1.data = tmp_st;
  
  an_st.data{iqid} = an_st_l1;
end % for iqid

an_st.meta.params = params;

end % function ensemble_enum_stats

%
% START OF VARIOUS SUB-FUNCTIONS
%

function report_stats_by_subject(data_st,params)
  col = set_var_col_const(data_st.vars);

end % report_stats_by_subject(an_st,params)

function report_stats_across_subjects(data_st,params)
  col = set_var_col_const(data_st.vars);

  % Deal with opening the file ID
  fid = ensemble_init_fid(params.tables);
  
  % Prepare variables for printing
  nsub = data_st.data{col.nsub};
  m = data_st.data{col.mean};
  sd = data_st.data{col.std};
  sem = sd/sqrt(nsub-1);
  
  if isfield(col,'mean_ttest')
    ttest_st = data_st.data{col.mean_ttest};
    ttest_cols = set_var_col_const(ttest_st.vars);
    pvalue = ttest_st.data{ttest_cols.p};
    tvalue = ttest_st.data{ttest_cols.stats}.tstat;
  else
    pvalue = NaN;
    tvalue = NaN;
  end
  
  if isfield(params,'question')
    qtxt = params.question.question_text;
    num_enum = length(params.question.enum_values);
    enum_str = sprintf('1=%s, %d=%s', ...
	params.question.enum_values{1}, ...
	num_enum, ...
	params.question.enum_values{num_enum});
    qtxt = sprintf('%s (%s):', qtxt, enum_str);
  else
    qtxt = '';
  end
    
  fprintf(fid,'%50s\tN\tMean\tSEM\tT\tprob\n','');
  fprintf(fid,'%50s\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.4f\n', qtxt, nsub, m, sem, tvalue, pvalue);
    
  if fid > 1
    fclose(fid);
  end
end % report_stats_across_subjects(an_st,params)

