function EEG = physio_filt_highlow(EEG,defs)

% high- and/or low-pass continuous data in an EEGlab struct using eegfilt
% 
%   EEG = physio_filt_highlow(EEG,defs)
% 
% REQUIRES
%   EEG - an EEGlab data structure
%   defs.physio_filt_highlow.high_cutoff - high pass filter cutoff in Hz
%   defs.physio_filt_highlow.low_cutoff - low pass filter cutoff in Hz
%   defs.physio_filt_highlow.filt_order
%   defs.physio_filt_highlow.channel - name of channel
% 
% RETURN
%   EEG - an EEGlab data structure, with the given signal having been
%   filtered
% 
% FB 2009.04.20

% try
  channel = defs.physio_filt_highlow.channel;
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
% catch
%   error('no target channel defined in defs.physio...artifact.channel\n');
% end
try lowcut = defs.physio_filt_highlow.low_cutoff; catch lowcut = []; end
try highcut = defs.physio_filt_highlow.high_cutoff; catch highcut = []; end
try fo = defs.physio_filt_highlow.filt_order; catch fo = []; end

if (~isempty(lowcut) || ~isempty(highcut))

  eeglab('initpaths');

  % Filter data  - do highpass and low-pass in separate stages to avoid
  % numerical problems
  if lowcut
    EEG.data(cidx,:) = eegfilt(EEG.data(cidx,:),EEG.srate,lowcut,0,0,fo);
  end
  if highcut
    EEG.data(cidx,:) = eegfilt(EEG.data(cidx,:),EEG.srate,0,highcut,0,[]);
  end
end
