function outdata = ensemble_export_cefa(indata,defs)

% outputs given dataset in a CEFA-friendly format, for factor analysis
% 
%   outdata = ensemble_export_cefa(indata,defs)
% 
% This function currently takes 'indata' (an N x M data matrix, N
% observations and M variables), calculates correlations between the
% variables, and writes the correlational data, as well as header
% information, to a text file in a format that CEFA (the Comprehensive
% Exploratory Factor Analysis program, Browne & Cudeck) can easily import.
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
%   - write out variable formats (raw, polychor corr)
%   - specify var names in 'defs'
%   - specify a confirmatory factor structure
% 
% REQUIRES
%   indata - N x M data matrix, N obs and M vars
%   defs.init_fid
%       defs.init_fid.write2file
%       defs.init_fid.print
%       defs.init_fid.fname
%       defs.init_fid.filemode
%       
% RETURNS
%   outdata - contents of the CEFA import file
% 
% FB 2010.03.30

outdata = [];

% init vars, init output file
fid = ensemble_init_fid(defs.init_fid);
nobs = size(indata,1);
nvar = size(indata,2);
data_type = 1; % 1=corr, 2=raw, 3=?, 4=?
var_names = cell(nvar,1);
for iv=1:nvar, var_names{iv} = sprintf('var%d',iv); end

% calculate correlation
r = corrcoef(indata);

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

% write out data
for j=1:nvar
  for k=1:j
    fprintf(fid,'%1.2f',r(j,k));
    if k<j, fprintf(fid,' '); end
  end
  fprintf(fid,'\n');
end

fclose(fid);
