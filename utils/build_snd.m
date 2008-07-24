function [waveform, T] = build_snd(freqs, amplitudes, phases, Fs, npts)

% function [waveform, T] = build_snd(freqs, amplitudes, phases, Fs, npts)
%
% Constructs a sound from the components specified in the freqs vector
% The freqs vector must be accompanied by a phases vector of equal length
% which specifies the starting phases of each frequency component in degrees

% Modification history:
% 12/95 	PJ 	Wrote function
% 1/15/96	PJ	Added amplitude vector
% 11/11/98      PJ      Fixed phase value
% 02/15/05      PJ      Phase value was incorrect.  Now really fixed

if max(freqs) > Fs/2
  error(['Maximum specified frequency <' int2str(max(freqs)) '> greater than nyquist (' int2str(Fs/2) ')']);
end

freqs = freqs(:)';
amplitudes = amplitudes(:)';
phases = phases(:)';

waveform = ones(npts,1)*amplitudes .* sin(2*pi*((0:npts-1)'/Fs*freqs + ones(npts,1)*phases/360));

T = (0:npts-1)/Fs;

if size(waveform,2) > 1
  waveform = sum(waveform,2);
end

return
