function [ses_mask, bad_chans, bad_trials] = bad_from_mask(ses_mask, bs)
%
% [ses_mask, bad_chans, bad_trials] = bad_from_mask(ses_mask, badspecs)
%
% Returns bad channels and modifies ses_mask based on current ses_mask and
% specs in badspecs
%   .prop_bad_chan  -- proportion of bad channels for trial to be declared bad
%   .prop_bad_trials  -- proportion of bad trials for channel to be declared bad
%
% ses_mask is a binary mask.  Ones indicate a good trial/channel.  Zero is
% bad. ses_mask is an M X N matrix in which the M rows represent trials and N
% columns represent channels.

% Modification history:
%
%  01/24/00 PJ
%  07/15/01 PJ Added refchan handling

if bs.add_refchan2mask
  ses_mask(:,end+1) = 1;
end

ses_mask = ~ses_mask;

bad_chans = find(sum(ses_mask) >= bs.prop_bad_trials*size(ses_mask,1));
bad_trials = find(sum(ses_mask,2) >= bs.prop_bad_chan*size(ses_mask,2));

ses_mask(bad_trials,:) = 1;
ses_mask(:,bad_chans) = 1;

ses_mask = ~ses_mask;

return

