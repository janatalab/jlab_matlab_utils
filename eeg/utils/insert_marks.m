function insert_marks(EEG)
% This function can be called to load a set of region marks from a mat file
% to the .winrej field of 'UserData' in the eegplot figure showing data
% that is being edited.
%
% insert_marks(EEG);
%

% 27Jun2012 Petr Janata

% Get the UserData from the current figure and make sure that it has a
% .winrej field

if nargin < 1
	error('Please provide EEG info structure')
end

g = get(gcf,'UserData');
if ~isfield(g, 'winrej')
	error(['No winrej field found in UserData.' ...
		'Make sure that the EEG data are showing in a window'])
end

% Specify a file to read rejection information from
[fpath,fstub,fext] = fileparts(EEG.filename);
outfname = fullfile(EEG.filepath, sprintf('%s_marks.mat',fstub));

% Make sure the file exists
if ~exist(outfname,'file')
	error('File with marks does not exist: %s\n', outfname)
end

fprintf('Loading marks from file: %s\n', outfname);
load(outfname, 'winrej');
fprintf('Found %d marks\n', size(winrej,1));

g.winrej = winrej;

% Write the information out to the UserData associated with the figure
% object
set(gcf,'UserData',g);
