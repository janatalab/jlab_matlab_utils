function outdata = ensemble_export_cefa(indata,defs)

% outputs given dataset in a CEFA-friendly format, for factor analysis
% 
%   outdata = ensemble_export_cefa(indata,defs)
% 
% This function currently takes 'indata' (in the format of either: an N x M
% data matrix, N observations and M variables; or an ensemble data struct),
% optionally calculates correlations between the variables, and writes
% either the correlational data, the raw data, or both as well as header
% information, to a text file in a format that CEFA (the Comprehensive
% Exploratory Factor Analysis program, Browne & Cudeck) can easily import.
% 
% NOTE: if your indata is in ensemble data struct format, please use
% defs.var_idxs to identify those numerical vars that you wish to output
% (to exclude the non-numerical vars). If indata format is ensemble data
% struct, it will skip all rows that contain non-numerical characters in a
% column value.
% 
% % % % % CEFA FORMAT
% 
% nobs nvars
% data_format
% 0
% 
% var_names
% (var names)
% 
% confirmatory_structure
% (confirmatory structure)
% 
% covariance_structure
% (covariance structure)
% 
% data
% 
% % % % % END CEFA FORMAT
% 
% FEATURE REQUESTS:
%   - accept ensemble data structures and other formats
%   - write out variable formats (polychor)
%   - specify var names in 'defs'
%   - specify a confirmatory factor structure
%   - deal with non-num values in ensemble data struct indata formats
%       (maybe automatically identify non-num columns?)
% 
% REQUIRES
%   indata - either an N x M data matrix (N obs and M vars) or an ensemble
%       data struct
%   defs.init_fid
%       defs.init_fid.write2file
%       defs.init_fid.print
%       defs.init_fid.fname
%       defs.init_fid.filemode
%   defs.var_idxs - the numerical indices of columns (or vars) of indata
%       that you would like to export to CEFA format. Default:
%       1:size(indata,2) or 1:length(indata.vars), depending on indata
%       format
%   defs.rev_score_idxs (optional) - a logical vector that indicates which
%       vars are negatively worded items. If this parameter is specified,
%       then the given variables will be reverse-scored. See:
%       defs.rev_score_max. NOTE: rev_score_idxs refer to the variables
%       that are present AFTER indata are masked by defs.var_idxs.
%       THEREFORE, if your indata has 100 variables, and you specify
%       def.var_idxs = [10 20 30 40 50], to specify the fact that the
%       original variables 20 and 50 need to be reverse scored, you would
%       set defs.rev_score_idxs = [2 5], NOT [20 50].
%   defs.rev_score_max (optional) - if defs.rev_score_idxs is specified,
%       then for each index in that parameter, the matching value in
%       rev_score_max, +1, will be subtracted from the given variable when
%       reversing that variable's score. If no rev_score_max is specified,
%       but rev_score_idxs is specified, then 1+max(var_values) will be
%       subtracted from the given values.
%   defs.parcel_idxs (optional) - cell array of vectors containing indices
%       of variables to include in each parcel. If parcel_idxs{1} = [1 3
%       5], then the first, third and fifth variables (after taking into
%       account defs.var_idxs) will be added together to construct parcel
%       1. if parcel_idxs{2} = [2 4 6] then the second, fourth, and sixth
%       variables (After taking into account defs.vars_idxs) will be added
%       together to construct parcel 2.
%   defs.var_names - cell array of strings identifying output variables. If
%       your input data format is ensemble data struct, and you do not
%       specify defs.var_names, it will be set to indata.vars{var_idxs}. If
%       your input data format is a numerical data matrix, var_names will
%       be set to sprintf('var%d',iv) for iv=var_idxs. If you specify
%       defs.parcel_idxs, var_names will be applied to  parcels.
%   defs.output_format {'raw','corr'} default: raw
%       
% RETURNS
%   outdata - contents of the CEFA import file
% 
% FB 2010.03.30
% FB 2010.05.25 - added support for raw output format (now the default),
%   both numerical matrix and ensemble data struct input, reverse scoring,
%   parcel creation, variable subsets, and variable naming

%% init vars
outdata = [];

% init output file
fid = ensemble_init_fid(defs.init_fid);

% get var_idxs
if isfield(defs,'var_idxs'), var_idxs = defs.var_idxs; end
if isfield(defs,'var_names'), var_names = defs.var_names; end

% get data
if isstruct(indata) && isfield(indata,'data')
  if ~exist('var_idxs','var')
    var_idxs = 1:length(indata.vars);
  end
  data = [indata.data{var_idxs}];
  if ~exist('var_names','var')
    var_names = {indata.vars{var_idxs}};
  end
elseif isnumeric(indata)
  if ~exist('var_idxs','var')
    var_idxs = 1:size(indata,2);
  end
  data = indata(:,var_idxs);
  if ~exist('var_names','var')
    var_names = cell(1,length(var_idxs));
    for iv=1:length(var_idxs)
      var_names{iv} = sprint('var%d',var_idxs(iv));
    end
  end
else
  error('unknown indata format\n');
end

% describe data
nobs = size(data,1);
nvar = size(data,2);
if isfield(defs,'output_format')
  fmt = defs.output_format;
else
  fmt = 'raw';
end

% set data_type, 1=corr, 2=raw, 3=?, 4=?
switch fmt
  case 'corr'
    data_type = 1;
  case 'raw'
    data_type = 2;
  otherwise
    error('unknown output format: %s\n',fmt);
end

%% reverse score?
if isfield(defs,'rev_score_idxs')
  rsis = defs.rev_score_idxs;
  nrsi = length(rsis);
  if isfield(defs,'rev_score_max')
    rsms = defs.rev_score_max;
  else
    rsms = [];
  end
  nrsm = length(rsms);
  for rsi=1:nrsi
    ridx = rsis(rsi);
    if rsi <= nrsm
      data(:,ridx) = (ones(nobs,1)*(rsms(rsi)+1)) - data(:,ridx);
    else
      data(:,ridx) = (ones(nobs,1)*(max(data(:,ridx))+1)) - data(:,ridx);
    end
  end
end % if isfield(defs,'rev_score_idxs

%% create parcels?
if isfield(defs,'parcel_idxs')
  npidx = length(defs.parcel_idxs);
  parcels = zeros(nobs,npidx);
  for k=1:npidx
    pidx = defs.parcel_idxs{k};
    parcels(:,k) = sum(data(:,pidx),2);
  end % for k=1:npidx
  data = parcels;
  nvar = npidx;
end % isfield(defs,'parcel_idxs

%% output headers

% write header info to file
fprintf(fid,'%d %d\n%d\n0\n\n',nobs,nvar,data_type);

% write variable names
if isempty(var_names)
  fprintf(fid,'0\n\n');
else
  fprintf(fid,'1\n');
  fprintf(fid,'%s\n\n',cell2str(var_names,' '));
end

% confirmatory factor structure
% % % currently not supported
fprintf(fid,'0\n0\n\n')

%% write out data

% data format?
switch fmt
  case 'raw'
    % write out data
    for k=1:nobs
      fprintf(fid,'%f ',data(k,:));
      fprintf(fid,'\n');
    end
  case 'corr'
    % calculate correlation
    r = corrcoef(data);

    % write out data
    for j=1:nvar
      for k=1:j
        fprintf(fid,'%1.2f',r(j,k));
        if k<j, fprintf(fid,' '); end
      end
      fprintf(fid,'\n');
    end
end


fclose(fid);

%% save outdata
outdata = ensemble_init_data_struct();
outdata.type = 'cefa_out';
if isfield(defs,'outDataName')
  outdata.name = defs.outDataName;
else
  outdata.name = outdata.type;
end
outdata.vars = var_names;
for k=1:length(outdata.vars)
  outdata.data{k} = data(:,k);
end

fprintf(1,'ensemble_export_cefa: DONE!\n');

