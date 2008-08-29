function [ses_mask] = artifact_edit(fname, fhdr, chdr, edit_file_type)
%  ses_mask=artifact_edit(fname, fhdr, chdr)
%  
%  Performs artifact editting of choice, according to the
%  proper edit_file_type.  Current options are:
%  		0)  Edit codes files as written by MacAverager
%		1)	Ramesh's MYSTERY format
%		n) 	YOUR format

%  Modification history 
%		Started July 7, 1995 by Ramesh Srinivasan and Petr Janata
%
%  MYSTERIES:see bad_channels_and_bad_trials in matlab directory

disp(['Starting edit code reading']);
ses_hdr_offsets_v; ave_hdr_offsets_v;

% Edit Code File Types
EDC = 0; MYSTERY = 1;
if nargin < 4
disp('Possible edit code file types are: ');
disp('  0) Edit codes files as written by MacAverager')
disp('  1) MYSTERY format')
edit_file_type = input('Select one: ')
end % if nargin < 4

if edit_file_type == EDC
 ses_mask = read_edc([fname '.edc'], fhdr, chdr);

elseif edit_file_type == MYSTERY

mys_fname = [fname '.mys'];
mys_fid = fopen(mys_fname, 'r');
ses_mask = zeros(sum(chdr(:,NTrials)), fhdr(NChan));
ses_mask = fread(mys_fid, [sum(chdr(:,NTrials)), fhdr(NChan)], 'int8');





end


ses_mask(:, fhdr(NChan) + 1) = ones(sum(chdr(:, NObs)), 1);

%replace average reference column in those rows that contain all
%zeros except for the avg. ref. column with zeros
%
replace_indices = find(sum(ses_mask') == 1)';
ses_mask(replace_indices, fhdr(NChan) + 1) = zeros(length(replace_indices),1);

