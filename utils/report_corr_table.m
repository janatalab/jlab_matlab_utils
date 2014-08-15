function report_corr_table(corr_st, params)
% Prints a correlation table
% INPUT:
% corr_st.r - matrix containing the correlations
% corr_st.p - matrix containing the corresponding p values
% corr_st.N - scalar indicating number of observations
%
% params.varLabels - cell array of variable names
% params.precision - number of significant digits to print
% params.paths.analyses - path to where the file should be saved
% params.fname - filename to associate with the file
% params.write2file - saved to file if true; otherwise reported to screen
% params.delim - string specifying the delimiter to use when writing the
% file

% 12Aug2014 Petr Janata

% Set any missing variables
if ~isfield(params,'delim')
  delim = '\t';
else
  delim = params.delim;
end

if ~isfield(params,'precision')
  precision = 2;
end

if ~isfield(params,'varLabels')
  varLabels = {};
else
  varLabels = params.varLabels;
end

% Perform sanity checks
nvars = size(corr_st.r,1);
nlabels = length(varLabels);
if nlabels ~= nvars
  error('Mismatch in number of variables (%d) and variable labels (%d)', nvars, nlabels)
end

% Get the file identifier
if strcmp(delim,',')
  fext = '.csv';
else
  fext = '.txt';
end
params.fname = fullfile(params.paths.analyses,[params.fstub fext]);
fid = ensemble_init_fid(params);

% Write the header line
fprintf(fid,delim);
for ivar = 1:nvars
  fprintf(fid,varLabels{ivar});
  if ivar < nvars
    fprintf(fid,delim);
  else
    fprintf(fid,'\n');
  end
end

% Write the data
fmt_str = sprintf('%%.%df%%s', precision);

for ivar = 1:nvars
  fprintf(fid,sprintf('%%s%s',delim),varLabels{ivar});
  for jvar = 1:nvars
    fprintf(fid,fmt_str,corr_st.r(ivar,jvar),prob2str(corr_st.p(ivar,jvar),0.05,'*',3,false));
    if jvar < nvars
      fprintf(fid,delim);
    else
      fprintf(fid,'\n');
    end
  end % for jvar
end % for ivar
