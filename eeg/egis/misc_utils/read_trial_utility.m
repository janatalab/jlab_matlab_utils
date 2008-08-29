function trialdata = read_trial_utility(fid, cell_data_offset, trial_num, nchan, npoints)
% trialdata = RD_ONETR_ALLCH(fileid, cell_data_offset, trial_num, nchan, npoints)
%
%   Reads the raw (untransformed) data from all channels for a single trial 
%   from a single cell in the EGIS session file referred to by the file_id.
%
%   Returns channels in columns

% Modification history:
%
%  2/15/95 PJ -- started work on module
%

trial_offset = cell_data_offset + (trial_num-1)*nchan*npoints*2;

fseek(fid, trial_offset, 'bof');

trialdata = fread(fid, [nchan, npoints],'int16');
trialdata = trialdata';
