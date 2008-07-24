function as = init_analysis_struct()
% Initializes an analysis structure for pulling data out of the mysql database
%
% as = init_analysis_struct()
%
% There are two basic sub-structures within an analysis structure. 
% 1) The .num structure contains numeric responses to questions. Most
%    importantly this includes enums.
% 2) The .txt structure contains text responses
%
% Both of these structures include additional fields with information for
% decoding the .data matrices.

% 08/19/05 Petr Janata

as.name = '';  % name of the analysis
as.forms = {};  % Comma-separated list of form names to aggregate into the analysis

as.num.subid = {}; % Vector of subject IDs
as.num.qid = [];  % Vector of question IDs
as.num.qtxt = {};  % Question text
as.num.vars = {};  % Output variable names
as.num.data = [];  % Output data matrix (nsubs X qid) for response data
as.num.datenum = [];  % Timestamps of response data
as.num.niter = []; % Number of iterations of each question
as.num.dfid = [];  % data format ID  - useful for enums
as.num.stimidx = []; % index into master stimulus list (nsubs x nreps)

as.txt.subid = {};
as.txt.qid = [];
as.txt.qtxt = {};
as.txt.vars = {};
as.txt.data = {};
as.num.datenum = [];  % Timestamps of response data
as.txt.niter = [];
as.txt.stimidx = [];

as.stims.ids = [];

as.params.stims = 1;
