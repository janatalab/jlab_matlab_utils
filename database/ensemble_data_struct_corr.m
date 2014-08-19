function out_st = ensemble_data_struct_corr(data_st,params)
% Calculates correlations between variables that are contained in the
% Ensemble data structure (data_st), and that are specified in
% params.corr.vars. 
%
% params.filt - filtering applied by ensemble_filter() prior to calculation
% of correlation matrix

% 11Mar2014 Petr Janata

% Initialize the output data structure
out_st = ensemble_init_data_struct;
out_st.type = 'correlation_matrix';
out_st.vars = {'corrvars','corrmtx','corrpval'};
ocols = set_var_col_const(out_st.vars);

% Make sure variables between which we want to calculate correlations are
% specified
if ~isfield(params,'corr') || ~isfield(params.corr,'vars') || ~iscell(params.corr.vars)
  error('Must specify variables to correlate as a cell array of strings in params.corr.vars')
else
  corrvars = params.corr.vars;
end

% Apply filtering params
if isfield(params,'filt')
  data_st = ensemble_filter(data_st,params.filt);
end

% Extract the data into a matrix
varmask = ismember(data_st.vars, corrvars);
extracted_vars = data_st.vars(varmask);
nvars = length(extracted_vars);
data = cell2mat(data_st.data(varmask));
N = size(data,1)-sum(any(isnan(data),2));

% Report any issues with NaNs
numnans = sum(isnan(data));
if any(numnans)
  fprintf('Encountered NaNs in the data:\n');
  for ivar = 1:length(extracted_vars);
    fprintf('\t%s: %d\n', extracted_vars{ivar}, numnans(ivar));
  end
end

% Calculate the correlation for complete rows
[corrdata, corrpval]= corrcoef(data,'rows','complete');
out_st.data{ocols.corrmtx} = corrdata;
out_st.data{ocols.corrpval} = corrpval;

out_st.data{ocols.corrvars} = extracted_vars;

% Make a plot if we want it
if isfield(params,'report') && isfield(params.report,'figs')
  if isfield(params.report.figs, 'plot')
    doPlot = params.report.figs.plot;
  else
    doPlot = true;
  end
  
  if isfield(params.report.figs, 'write2file')
    write2file = params.report.figs.write2file;
  else
    write2file = false;
  end
end

if doPlot
  figure
  load new_seismic
  imagesc(corrdata), axis xy image
  set(gca,'clim',[-1 1])
  colormap(new_seismic)
  colorbar
  
  set(gca,'ytick',1:nvars,'yticklabel', extracted_vars)
  
  set(gca,'xticklabel','')
  for ivar = 1:nvars
    text(ivar,0.4,strrep(extracted_vars{ivar},'_','\_'), ...
      'rotation',-90,'horizontalalign','left','verticalalign','middle');
  end
  
  title(sprintf('N = %d', N),'fontsize',18)
end

if write2file
  figfname = fullfile(params.paths.figures, [params.report.figs.fstub '.eps']);
  fprintf('Printing figure: %s\n', figfname);
  print(figfname, '-depsc')
end

return