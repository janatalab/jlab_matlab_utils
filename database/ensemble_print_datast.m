function outdata = ensemble_print_datast(indata,defs)

% prints the contents of a data struct to screen and/or to file
% 
%   outdata = ensemble_print_datast(indata,defs)
% 
% This simple script prints the contents of a given ensemble data structure
% to screen and/or to file. Given a single ensemble data structure, it will
% use "indata.vars" as the column names, and "indata.data" as the rows, and
% will print first the column names, then the rows, delimited by commas
% (default), or some other character (as defined by "defs").
% 
% Given nested ensemble data structures, or a cell array of data
% structures, it will iterate through each data structure and print the
% contents of each data structure. In this case, the function will recurse
% upon itself.
% 
% NOTE: if there are semi-nested data structures, such that in a given data
% structure there are some columns that are not nested data structures, and
% some that are, this funciton will not output the data in that particular
% structure.
% 
% NOTE: this function doesn't know how to deal with nested cells within
% indata.data columns. Currently, these are output as "MxN cell"
% 
% REQUIRES
%   indata.vars
%   indata.data
%   indata.name - if present, this will be used to identify data structure
%       output, especially in the case of nested data structures. if no
%       .name is found, but .type is found, .type will be used.
%   indata.parent_name - if the function is given nested data structures,
%       upon recursion, it will add this name to the resulting indata
%       structure, and the contents of this name will be used to identify
%       output at the given sub-struct level. indata.parent_name will be
%       defined thus:
% 
%       sprintf('%s.%s',parent_indata_name,child_indata_name)
% 
%   defs.delim (default: comma) - the character or set of characters that
%       will be used to delimit the output values
%   defs.fid - if nested data structures are provided as input, the file id
%       from the first iteration will be saved, passed on to recursive
%       calls to this function, and used for output over defs.export.
%   defs.export - parameters for ensemble_init_fid. if present, they will
%       be used to try to print data to file.
%       .print, .write2file, .fname, .filemode
%   defs.print2screen (default: 1) - print data to the screen (stdout)?
% 
% RETURNS
%   indata - in the form it was provided
% 
% FB <fbarrett@ucdavis.edu> 2010.10.18

%% set variables
outdata = indata;
ic = set_var_col_const(indata.vars);
nc = length(ic);

try delim = defs.delim; catch delim = ','; end
try print2screen = defs.print2screen; catch print2screen = 1; end

name = '';
if isfield(indata,'name')
  name = indata.name;
elseif isfield(indata,'type')
  name = indata.type;
elseif isfield(defs,'name')
  name = defs.name;
end
if ~isempty(name) && isfield(indata,'parent_name')
  name = sprintf('%s.%s',indata.parent_name,name);
end

if isfield(defs,'fid')
  fid = defs.fid;
  parent = 0;
elseif isfield(defs,'export')
  fid = ensemble_init_fid(defs.export);
  defs.fid = fid;
  parent = 1;
else
  fid = 0;
end

if ~fid && ~print2screen, error('no output specified'); end

%% check to see if there are nested data structures
nested = zeros(1,nc);
for k=1:nc
  if isstruct(indata.data{k})
    if isfield(indata.data{k},'data') && isfield(indata.data{k},'vars')
      % nested data structure is found
      nested(k) = 1;
      
      % recurse
      recdata = indata.data{k};
      if ~isempty(name), recdata.parent_name = name; end
      ensemble_print_datast(recdata,defs);
    else
      nested(k) = -1;
      msg = 'bad input structure found';
      if ~isempty(name), msg = sprintf('%s in %s',msg,name); end
      msg = sprintf('%s, SKIPPING',msg);
      warning(msg);
    end
  end % end if isstruct(indata.data{k
end % for k=1:nc

if any(nested) && ~all(nested)
  msg = 'some, but not all, nested data structures found';
  if ~isempty(name), msg = sprintf('%s in %s',msg,name); end
  msg = sprintf('%s, SKIPPING non-nested structures',msg);
  warning(msg);
  return
end

%% non-nested data structure, so print it out

% generate output format string
fmt = cell(nc,1);
for k=1:nc
  if isnumeric(indata.data{k})
    fmt{k} = '\d';
  else
    fmt{k} = '\s';
  end
end % for k=1:nc
fmtstr = sprintf('%s\n',cell2str(fmt,delim));

% print header lines
if ~isempty(name)
  namestr = sprintf('----------\nSource: %s\n----------\n\n',name); end
  if print2screen, fprintf(1,namestr); end
  if fid, fprintf(fid,namestr); end
end

varstr = sprintf('%s\n',cell2str(indata.vars,delim));
if print2screen, fprintf(1,varstr); end
if fid, fprintf(fid,varstr); end

% iterate over rows, print data
for k=1:length(indata.data{1})
  % generate data string for the given row
  input = cell(nc,1);
  for l=1:nc
    data = indata.data{l}(k);
    if isnumeric
      input = data;
    elseif iscell(data)
      if iscell(data{1})
        input = sprintf('%dx%d cell',size(data,1),size(data,2));
      else
        input = data{1};
      end
    end
  end
  inputstr = cell2str(input,delim);

  % print to screen and/or file
  if print2screen, fprintf(1,fmtstr,inputstr); end
  if fid, fprintf(fid,fmtstr,inputstr); end
end % for k=1:length(indata.data{1

%% clean up?
% print a few extra line endings
if print2screen, fprintf(1,'\n\n'); end
if fid, fprintf(1,'\n\n'); end

% if this call opened a file id, we should close it
if parent, fclose(fid); end
