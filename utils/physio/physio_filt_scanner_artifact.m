function EEG = physio_filt_scanner_artifact(EEG,defs)

% filters scanner artifact for a given physio channel
% 
%   EEG = physio_filt_scanner_artifact(EEG,defs)
% 
% This function uses EEGlab to epoch a given signal by TR, and then
% constructs a mean TR waveform and regresses the given signal by a design
% matrix containing a mean TR waveform regressor for each TR. It then
% removes the residuals from this regression from the signal of interest.
% 
% REQUIRES
%   EEG - an EEGlab data structure
%   defs.physio_filt_scanner_artifact.channel - name of channel to be
%       filtered, as it appears in EEG.chanlocs(:).labels
%   defs.TR - length of TR, in seconds
% 
% RETURNS
%   EEG - an EEGlab data structure, with the given signal having been
%   filtered
% 
% FB 2009.04.20

try TR = defs.TR; catch error('no TR defined in defs\n'); end
try channel = defs.physio_filt_scanner_artifact.channel; catch
    error('no target channel defined in defs.physio...artifact.channel\n');
end
if iscell(channel), channel = channel{1}; end

eeglab('initpaths');

cidx = strmatch(channel,{EEG.chanlocs(:).labels});
if isempty(cidx)
    error('target channel %s can not be found within EEG.chanlocs\n',...
        channel);
elseif length(cidx) > 1
    warning('many channels match target channel %s, using first one\n',...
        channel);
    cidx = cidx(1);
end

% create TR event channel
fprintf(1,'creating TR event channel\n');
tr_idxs = 1:TR*EEG.srate:EEG.pnts; % in samples
ntr = length(tr_idxs);
tr_events = [ones(ntr,1) (tr_idxs/EEG.srate)']; % in seconds
EEGe = pop_importevent(EEG,'append','no','event',tr_events,...
  'fields',{'type','latency'},'timeunit',1);

% extract TR epochs
fprintf(1,'extracting epochs\n');
EEGe = pop_epoch(EEGe,{},[0 TR],'epochinfo','yes');

% calculate mean TR waveform for each channel
fprintf(1,'calculating mean TR waveform for %s channel, ',channel);
W = [];
for it=1:EEGe.trials
  W = [W; EEGe.data(cidx,:,it)];
end
mW = mean(W); % mean waveform across TRs

% create regression design matrix
dM = zeros(EEG.pnts,EEGe.trials);
fprintf(1,'dM, ');
for it=1:EEGe.trials
  start = (EEGe.pnts*(it-1)+1);
  stop  = EEGe.pnts*it;
  dM(start:stop,it) = mW';
end

% are there any samples not accounted for at the end of dM?
if stop < size(dM,1)
  if (size(dM,1) - (stop - 1)) > length(mW)
    % left-over size is greater than the length of a TR
    error('there are less trials than there should be\n');
  else
    % fill out extra time at the end of the design matrix
    start = stop + 1;
    stop = size(dM,1);
    dM(start:stop,end+1) = mW(1:stop-start+1);
  end
end
dM(:,end+1) = 1; % must place a constant as the last regressor

fprintf(1,'regressing, ');
[b,bint,r,rint,stats] = regress(EEG.data(cidx,:)',dM);
EEG.data(cidx,:) = EEG.data(cidx,:) - r';
