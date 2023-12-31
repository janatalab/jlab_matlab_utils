function qtxt = mysql_get_qtxt(qids,varargin)
% Retrieves text associated with questions specified by the vector of qids
%
% qtxt = get_qtxt(qids,varargin)
%
% A flexible number of options can be passed into the function as tag/value
% pairs.
%
% 'conn_id', conn_id - the ID of the mysql connection - REQUIRED

% Parse the input arguments
narg = length(varargin);
for iarg = 1:2:narg
  switch varargin{iarg}
    case 'conn_id'
      conn_id = varargin{iarg+1};
    otherwise
      fprintf('get_qtxt: Unknown input argument: %s\n', varargin{iarg});
  end
end

% Check for valid connection to database
if ~exist('conn_id','var') || isempty(conn_id) || mysql(conn_id,'status')
  error('%s: Do not have a valid connection ID', mfilename);
end

qid_str = sprintf('%d,', qids);
qid_str = sprintf('(%s)', qid_str(1:end-1));

mysql_str = sprintf(['SELECT question_text FROM question ' ...
      'WHERE question_id IN %s;'], qid_str);

qtxt = mysql(conn_id,mysql_str);

return
