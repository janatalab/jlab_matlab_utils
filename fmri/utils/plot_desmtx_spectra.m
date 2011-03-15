function [data_st] = plot_spm_desmtx_spectra(flist,params)
% data_st = plot_desmtx_spectra(DesMtx,params)
%
% flist - List of SPM.mat files whose design matrices we want to plot
% spectra for.
%
% Parameters can be passed in for any spectral calculation method via a
% field in the params structure, where that field is named for the method,
% e.g. params.pwelch


% Make sure flist is a cellstr
if ~iscell(flist)
  flist = cellstr(flist);
end

% Loop over the SPM.mat files
nfiles = length(flist);
for ifile = 1:nfiles
  % Clear out any existing SPM structures
  clear SPM
  
  % Load the file
  fprintf('Working on file: %s\n', flist{ifile});
  load(flist{ifile})
  
  % Get the design matrix
  X = SPM.xX.X;
  
  % Get the column names
  col_names = SPM.xX.name;
  
  % Figure out which columns we are using
  try 
    use_cols = ismember(col_names, params.use_regress_names);
  catch
    use_cols = 1:length(col_names);
  end
      
  % calculate the power spectrum
  nrows = size(X,1);
  P = abs(fft(X(:,use_cols))).^2;
    
  % normalize by total power
  P = P./repmat(sum(P),size(P,1),1);
    
  data_st = P;
end % for ifile=

