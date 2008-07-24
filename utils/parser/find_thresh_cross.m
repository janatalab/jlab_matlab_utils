function [onsets, offsets] = find_thresh_cross(data,params)
% [onsets, offsets] = find_thresh_cross(data,params);
%
% Returns vectors containing the sample indices at which
% threshold values in the data vector were crossed. Will return both onsets and
% offsets. 
%
% params is a structure that governs how the threshold crossings are found
%
% .thresh - value that must be matched or exceeded
% .dir - direction to look at. If 'pos', then the algorithm look for values
%        that are more positive than the threshold value. If it is 'neg' then
%        the algorithm looks for values that are more negative than the
%        threshold value.
% .remove_mean - if true, the mean is removed from the data before the
%        thresholds are found. Note that the threshold value refers to the
%        zero-mean data.

% 07/20/06 Petr Janata

% Check for input arguments
error(nargchk(2,2,nargin))

onsets = [];
offsets = [];

% Make sure we have the required fields in the parameter structure
if ~isfield(params,'dir'), myfprintf('Missing .dir field\n'), return, end
if ~isfield(params,'thresh'), myfprintf('Missing .thresh field\n'), return, end

% Check to make sure that we are dealing with a vector
dims = size(data);
if all(dims > 1), myfprintf('Input must be a vector'), return, end

data = data(:);

% Remove the mean if desired
if isfield(params,'remove_mean') && params.remove_offset
  data = data-mean(data);
end

switch params.dir
  case {'pos','positive'}
    mask = data >= params.thresh;
  case {'neg','negative'}
    mask = data <= params.thresh;
  otherwise
    myfprintf(sprintf('Unknown direction: %s', params.dir));
    return
end

diff_vect = diff(mask);
onsets = find(diff_vect == 1)+1;
offsets = find(diff_vect == -1)+1;

function myfprintf(msg)
  fprintf('find_thresh_cross: %s\n', msg);
  return
  
