function save_marks(EEG)
% This function can be called to save the current set of region marks
% (stored in the .winrej field of 'UserData' in the eegplot figure) to a
% mat-file located in the same directory as the set file that is being
% edited 
%
% save_marks(EEG);
%

% 27Jun2012 Petr Janata

% Get the UserData from the current figure and make sure that it has a
% .winrej field

if nargin < 1
	error('Please provide EEG info structure')
end

g = get(gcf,'UserData');
if ~isfield(g, 'winrej')
	error('No winrej field found in UserData')
end

% Check to see if there are any marked regions to save
if isempty(g.winrej)
	fprintf('No marked regions to save\n');
	return
end

% Specify a file to write rejection information out to
[fpath,fstub,fext] = fileparts(EEG.filename);
outfname = fullfile(EEG.filepath, sprintf('%s_marks.mat',fstub));

winrej = g.winrej;

fprintf('Saving %d marks to file: %s\n', size(winrej,1), outfname);
save(outfname, 'winrej');

