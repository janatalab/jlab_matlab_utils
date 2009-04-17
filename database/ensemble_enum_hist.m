function an_st = ensemble_enum_hist(data_st,params)
% Plots the distributions of the response tallies in the categories of an enum.
%
% outdata = ensemble_enum_hist(data_st,params);
%
% Calculates and plots the distributions of the response tallies in the
% different categories of an enum. The function will operate on all qids that
% it finds in the input data that are of type enum
%
% Note: Currently, the script will not treat the same question appearing on
% different forms as a different instance of the question. If you don't want
% answers to the same question on different forms combined, you must filter the
% data to only process forms with unique question IDs.  This behavior may
% change in future versions.

% 01/28/07 Petr Janata - adapted from quest_resp_dist.m and an earlier start on
% ensemble_enum_hist which was much more oriented towards returning function
% handles to sub-functions that performed the various steps.
% 
% 02/17/07 PJ - adapted to handle either sessions, subjects, or stimuli as item
%               to be analyzed.
% 11/13/2007 FB - added error checking in plot_hist, to support blank enum 
% labels, such as in the case of data_format_id 142, 144 & 145
% 04/17/09 PJ - added output of plotted data to comma-separated files.

an_st = {};
na = 0;

try conn_id = params.ensemble.conn_id; catch conn_id = []; end

% Deal with some ambiguity in how reporting structures are defined
if isfield(params,'display') && isfield(params,'report')
  fprintf('%s: both display and report fields specified. Do not know which to use...\n', mfilename);
end
report_string_types = {'display','report'};
repstr = 'report';
for itype = 1:length(report_string_types)
  if isfield(params,report_string_types{itype})
    repstr = report_string_types{itype};
    break
  end
end

try do_plot = params.(repstr).figs.plot; catch do_plot = 0; end

% Apply any specified filtering to the input data
if isfield(params,'filt')
  fprintf('Applying filtering criteria\n')
  data_st = ensemble_filter(data_st, params.filt);
end

% Set the column constants
incol = set_var_col_const(data_st.vars);

%
% Gather metadata on the questions
%

% Figure out whether questions are referred to by their qids or composite qids
if isfield(incol,'question_id')
  qids = unique(data_st.data{incol.question_id});
elseif isfield(incol,'compqid')
  qids = fix(unique(data_st.data{incol.compqid}));
else
  fprintf('Did not find question_id column in the input data\n');
  return
end

qinfo = mysql_extract_metadata('table','question', ...
  'question_id',qids, ...
  'conn_id', conn_id);

% Generate the compqid for each row in the data if necessary
compqid_vect = [];
if isfield(incol,'compqid')
  compqid_vect = data_st.data{incol.compqid};
elseif isfield(incol,'subquestion')
  compqid_vect = make_compqid(data_st.data{incol.question_id}, ...
      data_st.data{incol.subquestion});
end

% Figure out which of the questions in the qinfo structure are enums
qinfo_enum_mask = ismember({qinfo.type},'enum');
qinfo_enum_idxs = find(qinfo_enum_mask);
nqid = length(qinfo_enum_idxs);

% Figure out which of the questions are bitmasks that allow for selection of
% multiple values (checkbox as opposed to radiogroup)
qinfo_bitmask_mask = ismember({qinfo.html_field_type},'checkbox');

qinfo_notbitmask_mask = qinfo_enum_mask & ~qinfo_bitmask_mask;

% Create masks for all of the response data
enum_compqids = [qinfo(qinfo_enum_mask).compqid];
enum_mask = ismember(compqid_vect, enum_compqids);
  
notbitmask_compqids = [qinfo(qinfo_notbitmask_mask).compqid];
notbitmask_mask = ismember(compqid_vect, notbitmask_compqids);

% Copy all of the enum data to the data_vect
data_vect = zeros(size(compqid_vect));
data_vect(enum_mask) = data_st.data{incol.response_enum}(enum_mask);

% Convert the data for all of the non-bitmask enums to category indices
data_vect(notbitmask_mask) = enum2data(data_vect(notbitmask_mask));

% Figure out which variable is being used as the item for the item
if isfield(params,'item_var')
  item_str = params.item_var;
else
  item_str = 'session_id';
end

switch item_str
  case 'session_id'
    item_type = 'session';
  case 'subject_id'
    item_type = 'subject';
  case 'stimulus_id'
    item_type = 'stimulus';
  otherwise
    fprintf('%s: Unknown item variable: %s\n', mfilename, item_str);
end

% Precalculate the masks for the item variable we are using
[item_mask_mtx, itemids] = make_mask_mtx(data_st.data{incol.(item_str)});
nitems = length(itemids);

%
% Set stuff up for writing to a file, if that's what we're going to do.
%
fid = ensemble_init_fid(params.(repstr).tables);

%
% Loop over all of the unique question/subquestion combinations or compqids
%
for iqid = 1:nqid
  curr_qid_idx = qinfo_enum_idxs(iqid);
  
  % Copy the question info over to the display parameter structure in case we
  % are going to display some of the data
  params.(repstr).qinfo = qinfo(curr_qid_idx);

  % Make sure this question is an enum
  if ~enum_mask(iqid)
    continue
  end

  % Get the enum values
  enum_values = qinfo(curr_qid_idx).enum_values;
  ncat = length(enum_values);
  
  % Make a mask for the data corresponding to this question
  qid_mask = ismember(compqid_vect,qinfo(curr_qid_idx).compqid);
  
  %
  % Perform analyses that retain item identity
  %
  na = na+1;
  
  an_st{na} = init_analysis_struct;
  an_st{na}.type = sprintf('enum_category_by_%s', item_type);
  an_vars = {item_str,'count','prop','nresp'};
  an_st{na}.vars = an_vars;
  an_cols = set_var_col_const(an_vars);
  by_item_an_cols = an_cols;
  by_item_an_idx = na;
  nr = 0;  % set the tally of reports for this analysis to zero
  
  % Initialize output variables
  for ia = 1:length(an_vars)
    switch an_vars{ia}
      case {'count','prop'}
	an_st{na}.data{ia} = zeros(nitems,ncat);
      case item_str
	an_st{na}.data{ia} = itemids;
      case {'nresp'}
	an_st{na}.data{ia} = zeros(nitems,1);
    end
  end
  
  % Prepare an output data matrix that tallies the number of times each
  % category was selected for each item
  
  for iitem = 1:nitems
    item_mask = item_mask_mtx(:,iitem);

    % Tally the number of responses the made for this item and this question
    nresp = sum(item_mask&qid_mask);
    an_st{na}.data{an_cols.nresp}(iitem) = nresp;
    
    if ~nresp
      continue
    end
    
    switch qinfo(curr_qid_idx).html_field_type
      case {'radiogroup','menu'}
	an_st{na}.data{an_cols.count}(iitem,:) = hist(data_vect(item_mask&qid_mask),1:ncat)';
      otherwise
	bitmask = data2bitmask(data_vect(item_mask&qid_mask),ncat);
	if size(bitmask,1) > 1
	  sumvect = sum(bitmask);
	else
	  sumvect = bitmask;
	end
	an_st{na}.data{an_cols.count}(iitem,:) = sumvect;
    end
  end % for iitem

  nresp_mtx = repmat(an_st{na}.data{an_cols.nresp},1,ncat);
  nresp_mtx(nresp_mtx == 0) = NaN;
  an_st{na}.data{an_cols.prop} = an_st{na}.data{an_cols.count} ./ nresp_mtx;

  an_st{na}.meta.question = qinfo(curr_qid_idx);

  %
  % Perform any reports of the data in which reporting is done on a per-item basis
  %
  nr=nr+1;
  an_st{na}.report{nr}.type = 'distrib_by_cat';
  an_st{na}.report{nr}.figfun = @plot_distrib_by_cat;
  if do_plot
    an_st{na}.report{nr}.figs = ...
	plot_distrib_by_cat(an_st{na}.data{an_cols.prop}, params.(repstr)); 
  end
  
  fprintf(fid,'\nQuestion (%1.2f): %s\n', ...
      qinfo(curr_qid_idx).compqid, qinfo(curr_qid_idx).question_text);
  fprintf(fid,'%s\t%s\n', item_str, cell2str(enum_values,'\t'));
  for iitem = 1:nitems
    switch item_str
      case {'session_id','stimulus_id'}
	itemval_str = sprintf('%d',itemids(iitem));
      case 'subject_id'
	itemval_str = itemids{iitem};
      otherwise
	itemval_str = itemids{iitem};
    end
    
    enum_str = sprintf('\t%d', an_st{na}.data{an_cols.count}(iitem,:));
    fprintf(fid,'%s%s\n', itemval_str, enum_str);
  end
  
  %
  % Now perform the set of analyses that collapse across items
  %
  na = na+1;
  
  an_st{na} = init_analysis_struct;
  an_st{na}.type = sprintf('enum_category_across_%s', item_type);
  an_vars = {'prop','count'};
  an_st{na}.vars = an_vars;
  an_cols = set_var_col_const(an_vars);
  nr = 0;

  for ivar = 1:length(an_vars)
    an_type = an_vars{ivar};
    switch an_type
      case 'prop'
	data_col = by_item_an_cols.prop;
	plot_title = 'Average proportion of responses in each category';
      case 'count'
	data_col = by_item_an_cols.count;
	plot_title = 'Overall counts within in each category';
	
    end % switch an_type
    data = an_st{by_item_an_idx}.data{data_col};
    
    % Create a subanalysis structure so that we can accumulate distributional
    % statistics for each category
    sa_st = init_analysis_struct;
    sa_st.type = 'distrib_stats';
    sa_vars = {'mean','std','min','max','nitems'};
    sa_cols = set_var_col_const(sa_vars);
    sa_st.vars = sa_vars;
    
    % Loop over all of the things we want to calculate
    for jvar = 1:length(sa_vars)
      sa_str = sa_vars{jvar};
      switch sa_str
	case 'mean'
	  if any(isnan(data)) fun = @nanmean; else fun = @mean; end
	  result = fun(data);
	case 'nitems'
	  result = sum(~isnan(data(:,1)));
	case 'std'
	  if any(isnan(data)) fun = @nanstd; else fun = @std; end
	  result = fun(data);
	case 'min'
	  if any(isnan(data)) fun = @nanmin; else fun = @min; end
	  result = fun(data);
	case 'max'
	  if any(isnan(data)) fun = @nanmax; else fun = @max; end
	  result = fun(data);
	  
      end % switch sa_str
      
      % Assign the result
      sa_st.data{sa_cols.(sa_str)} = result;
    end % for jvar
    
    % Assign the sub-analyses to the main analysis
    an_st{na}.data{ivar} = sa_st;
    
    %
    % Perform any reports of the data in which reporting is done on an across-items basis
    %
    nr=nr+1;
    an_st{na}.report{nr}.type = 'distrib_by_cat';
    an_st{na}.report{nr}.figfun = @plot_hist;
    if do_plot
      switch an_type
	case 'prop'
	  plot_title = 'Average proportion of responses in each category';
	  ylim = [0 1];
	case 'count'
	  plot_title = 'Overall counts within in each category';
	  ylim = []; % do it based on max at time of display
      end
      params.(repstr).figs.ylim = ylim;
      params.(repstr).figs.title = plot_title;
      an_st{na}.report{nr}.figs = plot_hist(an_st{na}.data{an_cols.(an_type)},params.(repstr)); 
    end % if do_plot
  end % for ivar=
  an_st{na}.meta.question = qinfo(curr_qid_idx);
  
  fprintf(fid,'\nQuestion (%1.2f): %s\n', ...
      qinfo(curr_qid_idx).compqid, qinfo(curr_qid_idx).question_text);
  for ienum = 1:length(qinfo(curr_qid_idx).enum_values)
    fprintf(fid,'%30s:\t%1.2f\n', qinfo(curr_qid_idx).enum_values{ienum}, ...
	an_st{na}.data{an_cols.prop}.data{sa_cols.mean}(ienum));
  end
end % for iqid

if fid > 1
  fclose(fid);
end

end % function ensemble_enum_hist

function an_st = init_analysis_struct
  an_st.type = '';
  an_st.vars = {};
  an_st.data = {};
  an_st.meta = [];
  
end % function an_st = init_analysis_struct

function out_pp = plot_hist(data_st,params)

% FB 11/13/2007 - added error checking, to support blank enum labels, such
% as in the case of data_format_id 142, 144 & 145

  pp = params.figs;

  % Parse some of the input parameters
  try use_fig = params.use_fig;
  catch 
    try use_fig = pp.start_fignum; catch use_fig = []; end
  end
  
  try use_axes = params.use_axes; catch use_axes = []; end
  
  try labelfontsize = params.axislabelfontsize; catch labelfontsize = 9; end
  
  col = set_var_col_const(data_st.vars);
  
  ncat = length(params.qinfo.enum_values);
  nitems = data_st.data{col.nitems};

  max_chars_per_line = 50;
  
  %
  % Generate a figure that summarizes the distributions across categories
  %
  if isempty(use_fig)
    figure
  else
    figure(use_fig)
  end
  
  nax = 0;
  
  nax=nax+1;
  if ~isempty(use_axes)
    axes(use_axes)
  end
  
  % Check to see if we need to sort the data
  try plot_sorted = pp.sort; catch plot_sorted = 0; end
  
  if plot_sorted
    try sortdir = pp.sort_dir; catch sortdir = 'descend'; end
    [sorted_data, sorted_idxs] = sort(data_st.data{col.mean},sortdir);
  else
    sorted_idxs = 1:length(data_st.data{col.mean});
  end
  
  % Create output variables
  sorted_data = data_st.data{col.mean}(sorted_idxs);
  sorted_stderr = sorted_data/sqrt(data_st.data{col.nitems});
  
  % Plot the data
  bs = bar(1:ncat, sorted_data);
  add_errorbars(bs,sorted_stderr');
  set(gca,'xtick',[],'xlim',[0 ncat+1])
  set(gca,'activepositionproperty','outerposition')

  ax(nax) = gca;
  
  % Set axes formatting
  
  % Add xtick labels as text objects so that we can rotate them
  label_data = {};
  for icat = 1:ncat
    % FB 11/13/2007 added this conditional, since empty enum values create
    % an empty cell return from linewrap, which in turn is illegal in
    % cell2str, which finally kills any analysis in its tracks
    if(length(linewrap(params.qinfo.enum_values{sorted_idxs(icat)},18)))
      label_text = cell2str(linewrap(params.qinfo.enum_values{sorted_idxs(icat)},18),'\n');
    else
      label_text = ' ';
    end
    text(icat,0,label_text,'rotation',-90, ...
      'horizontalalign','left', ...
      'verticalalign','middle', ...
      'fontsize', labelfontsize);
    label_data{icat} = label_text;
  end
  
  try ylim = pp.ylim; catch ylim=[]; end
  if isempty(ylim)
    ylim = [0 1.2*max(data_st.data{col.mean}(sorted_idxs))];
  end
  set(gca,'ylim',ylim)
  
  try ylabel_str = pp.ylabel; catch ylabel_str = ''; end
  ylabel(ylabel_str)
  
  if ~isfield(pp,'title')
    pp.title = '';
  end
  title(pp.title)
  
  if ~isfield(pp,'add_nitems') || pp.add_nitems
    th = add_nitems_txt(data_st.data{col.nitems});
  end  

  % Add a page header
  try add_pagehdr = params.add_pagehdr; catch add_pagehdr = 1; end
  if add_pagehdr
    pp.pagehdr.title = sprintf(['%s\nQuestion (%1.2f): %s'], ...
	pp.title, params.qinfo.compqid, ...
	cell2str(linewrap(params.qinfo.question_text,max_chars_per_line),'\n'));
    nax = nax+1;
    ax(nax)=add_fighdr(pp.pagehdr);
    pp.fig.title = pp.pagehdr.title;
  end
  
  if isfield(pp,'write2file') && pp.write2file
    fprintf('Printing figure to file: %s\n', pp.figfname);
    print(pp.figfname, pp.printargs{:})
  end
  
  % Write the data that were plotted to a comma-delimited text file
  try data2file = pp.data2file; catch data2file = true; end
  if data2file && ~isempty(pp.figfname)
    [fpath,fname,fext] = fileparts(pp.figfname);
    datafname = fullfile(fpath,[fname '.csv']);
    % open the file for writing
    fid = fopen(datafname,'wt');
    if fid == -1
      fprintf('FAILED to open %s for writing\n', datafname);
    else
      fprintf('Writing figure data to: %s\n', datafname);
      fprintf(fid,'Category,Data,StdErr\n');
      for icat = 1:length(label_data)
        fprintf(fid,'%s,%1.6f,%1.6f\n', label_data{icat}, sorted_data(icat), sorted_stderr(icat));
      end
      fclose(fid);
    end
  end
  
  pp.fignum = gcf;
  pp.axes = ax;
  
  out_pp = pp;
end % plot_hist()

function pp = plot_distrib_by_cat(data_st,params)
  pp = params.figs;
  max_chars_per_line = 50;

  % Parse some of the input parameters
  try use_fig = params.use_fig;
  catch 
    try use_fig = pp.start_fignum; catch use_fig = []; end
  end
  
  try use_axes = params.use_axes; catch use_axes = []; end
  
  try labelfontsize = params.axislabelfontsize; catch labelfontsize = 9; end

  if ~isstruct(data_st)
    tmp.vars = {'prop'};
    tmp.data = {data_st};
    data_st = tmp;
  end
  col = set_var_col_const(data_st.vars);
  ncat = length(params.qinfo.enum_values);
  nitems = size(data_st.data{col.prop},1);

  %
  % Generate a figure that shows distributions of responses across items
  % separately for each enum category
  %
  if isempty(use_fig)
    figure
  else
    figure(use_fig)
    clf
  end

  nax = 0;
  ax = [];
  
  try ncol = params.ncol; catch ncol = 3; end
  try tick_interval = pp.tick_interval; catch tick_interval = 0.1; end
  hist_scale = 0:tick_interval:1;
  
  for icat = 1:ncat
    hist_vals = hist(data_st.data{col.prop}(:,icat),hist_scale);
    
    try ylim = params.ylim; catch ylim=[]; end
    if isempty(ylim)
      ylim = [0 1.2*max(hist_vals)];
    end
    
    nax = nax+1;
    ax(nax) = subplot(fix(ncat/ncol)+rem(ncat,ncol),ncol,icat);
    bar(hist_scale,hist_vals)
    set(gca,'xtick',0:tick_interval:1,'xlim',[-0.05 1.05])
    set(gca,'ylim',ylim)
    title_str = params.qinfo.enum_values{icat};
    if strcmpi(title_str,'parent')
      title_str = [title_str ' '];
    end
    title(title_str)
    if ~isfield(pp,'add_nitems') || pp.add_nitems
      th = add_nitems_txt(nitems);
    end  
  end % for icat=
  
  pp.fig.fignum = gcf;

  if ~isfield(pp,'title')
    pp.title = '';
  end

  pp.pagehdr.title = sprintf(['%s\nQuestion (%1.2f): %s'], pp.title, ...
      params.qinfo.compqid,...
      cell2str(linewrap(params.qinfo.question_text,max_chars_per_line),'\n'));
  nax = nax+1;
  ax(nax)=add_fighdr(pp.pagehdr);

  pp.fig.title = pp.pagehdr.title;
  pp.fig.axes = ax;

  if isfield(pp,'write2file') && pp.write2file
    fprintf('Printing figure to file: %s\n', pp.figfname);
    print(pp.figfname, pp.printargs{:})
  end

  out_pp = pp;
end % function pp = plot_distrib_by_cat(data_st,params)

function th = add_nitems_txt(nitems,params)
  th = text(0.95,0.95,sprintf('N=%d', nitems), ...
      'units','norm', ...
      'horizontalalign','right', ...
      'verticalalign', 'top');
end % function add_nitems_txt(nitems,params)
