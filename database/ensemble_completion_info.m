function [out_st] = ensemble_completion_info(data_st,params)
% Returns completion statistics for a given data set.
% [out_st] = ensemble_completion_info(data_st,params);
%
% Returns various completion statistics that can be garnered from the data
% table provided in data_st.  The fields in the params structure control the
% behavior of this function.
%
% .filt - applies filtering based on fields in the exclude and include structures

out_st.type = 'ensemble_completion_info';
out_st.vars = {};
out_st.data = {};

% Apply any specified filtering to the input data
if isfield(params,'filt')
  data_st = ensemble_filter(data_st, params.filt);
end

% Get the column constants
INCOL = set_form_col_const(data_st.vars);

% By default, organize information by session_id
if isempty(INCOL.SESS_ID)
  fprintf('ensemble_completion_info: Could not find session_id variable\n');
  return
end

% Generate a list of output variables. This should be done dynamically based on
% requested information to compute, but for now it is a fixed list
out_vars = {'session_id','start_time','stop_time','elapsed_time','form_id'};
outcol = set_var_col_const(out_vars);

out_st.vars = out_vars;
sesslist = unique(data_st.data{INCOL.SESS_ID});
out_st.data{outcol.session_id} = sesslist;

% Pull the data for each session
nsess = length(sesslist);

for isess = 1:nsess
  % Generate a session mask
  sess_mask = data_st.data{INCOL.SESS_ID} == sesslist(isess);
  sess_idxs = find(sess_mask);
  
  % Get the start time
  start_time = min(data_st.data{INCOL.DATE_TIME}(sess_mask));
  out_st.data{outcol.start_time}(isess) = start_time;

  % Get the stop time
  [stop_time, stop_idx] = max(data_st.data{INCOL.DATE_TIME}(sess_mask));
  out_st.data{outcol.stop_time}(isess) = stop_time;
  
  % Get the elapsed time
  elapsed_time = etime(datevec(stop_time),datevec(start_time));
  out_st.data{outcol.elapsed_time}(isess) = elapsed_time;
  
  % Get the form ID of the last form
  last_form_id = data_st.data{INCOL.FORM_ID}(sess_idxs(stop_idx));
  out_st.data{outcol.form_id}(isess) = last_form_id;

end % for isess=

% Calculate various summary statistics
report.min_time = min(out_st.data{outcol.elapsed_time});
report.max_time = max(out_st.data{outcol.elapsed_time});
report.mean_time = mean(out_st.data{outcol.elapsed_time});
report.std_time = std(out_st.data{outcol.elapsed_time});
report.median_time = median(out_st.data{outcol.elapsed_time});

%
% Generate plots and tables if desired
%
if isfield(params,'display')
  if isfield(params.display, 'figs')
    figtypes = fieldnames(params.display.figs);
    ntypes = length(figtypes);
    for itype = 1:ntypes
      curtype = figtypes{itype};
      st = params.display.figs.(curtype);
      switch curtype
	case 'hist'  % histogram of elapsed times
	  try hist_min = st.min_time; catch hist_min = 0; end
	  try hist_max = st.max_time; catch hist_max = report.max_time; end
	  try nbins = st.nbins; catch nbins=20; end
	  try doplot = st.plot; catch doplot=1; end
	  try tick_interval = st.tick_interval; catch tick_interval=15; end
	  
	  histscale = linspace(hist_min,hist_max,nbins);
	  histdata = histc(out_st.data{outcol.elapsed_time}, histscale);
	  histscale = histscale/60;
	  report.hist.xscale = histscale;
	  report.hist.data = histdata;
	  
	  if doplot
	    figure
	    bar(histscale,histdata)
	    set(gca,'xtick',min(histscale):tick_interval:max(histscale), ...
		'xlim',[min(histscale) max(histscale)])
	    report.hist.axes = gca;
	  end
      end % switch curtype
    end % for itype = 1:ntypes
  end % if isfield(params.display, 'figs')
  
  if isfield(params.display, 'tables')
    if isfield(params.display.tables,'group_stats') && ...
	  params.display.tables.group_stats
      fprintf('Completion time statistics\n');
      fprintf('Minimum elapsed time (min): %1.2f\n', report.min_time/60);
      fprintf('Maximum elapsed time (min): %1.2f\n', report.max_time/60);
      fprintf('Mean elapsed time(min): %1.2f\n', report.mean_time/60);
      fprintf('Std. dev elapsed time(min): %1.2f\n', report.std_time/60);
      fprintf('Median elapsed time(min): %1.2f\n', report.median_time/60);
      
    end
  end % if isfield(params.display, 'tables')
end % if isfield(params,'display')

out_st.report = report;
