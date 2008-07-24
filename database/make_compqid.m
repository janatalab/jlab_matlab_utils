function compqid = make_compqid(qids,subqids, num_places);
% Makes a vector of "composite question IDs" from given question IDs.
%
% compqid = make_compqid(qids,subqids)
%
% Makes a vector of "composite question IDs" from the vectors of question IDs
% and subquestion numbers.
%
% The format for the composite question ID is question.subquestion
%
% num_places specifies how many places following the decimal point should be
% used for representing the subquestion ID. When num_places=1, there is a limit
% on 9 subquestions. For num_places=2 there is a limit of 99, etc.  Currently,
% the Ensemble-wide assumption is that num_places = 2.
%

% 01/28/07 Petr Janata

if nargin < 3
  num_places = 2;
end

compqid = qids + subqids * 10^-num_places;

