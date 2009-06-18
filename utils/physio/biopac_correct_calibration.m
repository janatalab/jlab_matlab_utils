function EEG = biopac_correct_calibration(EEG,defs)

% correct a biopac waveform collected with incorrect calibration
% 
%   EEG = biopac_correct_calibration(EEG,defs)
% 
% This script currently only corrects scr that was not calibrated (i.e. the
% calibration settings were the default settings).
% 
% REQUIRES
%   EEG - an EEGlab data structure
%   defs.biopac_correct_calibration.gain - gain setting of amplifier
%   defs.biopac_correct_calibration.channel - name of channel
% 
% RETURN
%   EEG - an EEGlab data structure, with the given signal having been
%   filtered
% 
% FB 2009.06.10

channel = defs.biopac_correct_calibration.channel;
if iscell(channel), channel = channel{1}; end
cidx = strmatch(channel,{EEG.chanlocs(:).labels});
if isempty(cidx)
  error('target channel %s can not be found within EEG.chanlocs\n',...
      channel);
elseif length(cidx) > 1
  warning('many channels match target channel %s, using first one\n',...
      channel);
  cidx = cidx(1);
end

switch channel
  case {'scr','gsr'}
    try gain = defs.biopac_correct_calibration.gain; catch gain = []; end

    if ~isempty(gain)
      EEG.data(cidx,:) = EEG.data(cidx,:)*gain;
    else
      error('no gain provided for SCR calibration');
    end
  otherwise
    warning('unknown channel %s',channel);
end
