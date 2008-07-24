function fh = ensemble_plotfun(plot_name)
% Returns a function handle to one of the plotting sub-functions.
% fh = ensemble_plotfun(plot_name);
%
% Returns a function handle to one of the plotting sub-functions.  Calling the
% function without any arguments returns the list of available functions
%

% 02/18/07 Petr Janata

plotfun_names = {'plot_hist','plot_distrib_by_cat'};
fh = [];
if nargin < 1
  fprintf('Available plotting functions:\n');
  fprintf('\t%s\n', cell2str(plotfun_names,'\n'));
  return
end

if isempty(strmatch(plot_name,plotfun_names))
  fprintf('Unable to find desired function: %s\n', plot_name);
else
  fh = str2func(plot_name);
end

end % function fh = ensemble_plotfun(plot_name)

%
% Main plotting functions
%
function out_pp = plot_hist(data_st,params)
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
  
  % Plot the data
  bs = bar(1:ncat,data_st.data{col.mean});
  add_errorbars(bs,data_st.data{col.std}'/sqrt(data_st.data{col.nitems}));
  set(gca,'xtick',[],'xlim',[0 ncat+1])
  set(gca,'activepositionproperty','outerposition')

  ax(nax) = gca;
  
  % Set axes formatting
  
  % Add xtick labels as text objects so that we can rotate them
  for icat = 1:ncat
    label_text = cell2str(linewrap(params.qinfo.enum_values{icat},18),'\n');
    text(icat,0,label_text,'rotation',-90, ...
	'horizontalalign','left', ...
	'verticalalign','middle', ...
	'fontsize', labelfontsize);
  end
  
  try ylim = pp.ylim; catch ylim=[]; end
  if isempty(ylim)
    ylim = [0 1.2*max(data_st.data{col.mean})];
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
    pp.pagehdr.title = sprintf(['%s\nQuestion (%1.2f): %s'], pp.title, params.qinfo.compqid, ...
	cell2str(linewrap(params.qinfo.question_text,max_chars_per_line),'\n'));
    nax = nax+1;
    ax(nax)=add_fighdr(pp.pagehdr);
    pp.fig.title = pp.pagehdr.title;
  end
  
  if isfield(pp,'write2file') && pp.write2file
    fprintf('Printing figure to file: %s\n', pp.figfname);
    print(pp.figfname, pp.printargs{:})
  end
  
  pp.fignum = gcf;
  pp.axes = ax;
  
  out_pp = pp;
end % plot_hist()


function pp = plot_distrib_by_cat(data_st,params)
  pp = params.figs;
  max_chars_per_line = 50;

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
  if isfield(pp,'start_fignum')
    figure(pp.start_fignum+nfig)
  else
    figure
  end
  
  nax = 0;
  ax = [];
  
  try ncol = pp.ncol; catch ncol = 3; end
  try tick_interval = pp.tick_interval; catch tick_interval = 0.1; end
  hist_scale = 0:tick_interval:1;
  
  for icat = 1:ncat
    nax = nax+1;
    ax(nax) = subplot(fix(ncat/ncol)+rem(ncat,ncol),ncol,icat);
    hist_vals = hist(data_st.data{col.prop}(:,icat),hist_scale);
    bar(hist_scale,hist_vals)
    set(gca,'xtick',0:tick_interval:1,'xlim',[-0.05 1.05])
    set(gca,'ylim',[0 1.1*max(hist_vals)])
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

  pp.pagehdr.title = sprintf(['%s\nQuestion (%1.2f): %s'], pp.title, params.qinfo.compqid,...
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


%
% Internal support functions
%

function th = add_nitems_txt(nitems,params)
  th = text(0.95,0.95,sprintf('N=%d', nitems), ...
      'units','norm', ...
      'horizontalalign','right', ...
      'verticalalign', 'top');
end % function add_nitems_txt(nitems,params)
