function norm_dir(indir,instub,outdir,outstub)
% norm_dir(dirname,instub,outstub);
%
% Performs global intensity normalization
%
% dirname -- directory containing .img files to be normalized
% instub -- root of filenames to be included in normalization,
%           e.g. 01aprFOOL_s02
% outdir -- directory to dump normalized files into
% outstub -- root name of output files
%

% 05/18/00 PJ

if nargin < 2
  error('Too few args to norm_dir')
end

if nargin < 3
  outdir = indir;
  outstub = sprintf('gn%s', instub);
elseif nargin < 4
  outstub = sprintf('gn%s', instub);
end

V = load_dir(indir,instub);

nV = norm_intensity(V);

write_dir(nV, outdir, outstub);

