function [scaled] = scale_wvf(wvf,dB,nbits,is_wav)
% [scaled] = scale_wvf(wvf,dB,nbits,is_wav);
% 
% Scales a waveform to a desired decibel level given a resolution of nbits.
% Default is 16 bits.
% Removes DC component
% If is_wav is set to true, the result is scaled on the range of -1 to 1 to be
% compatible with .WAV format [default=false]
%
% If wvf contains multiple columns, each column is scaled independently.

% 01/05/05 PJ
% 06/21/06 PJ Fixed a problem in that the scaling factor was being calculated
%             based on the desired peak value and the rms of the waveform.

% are we dealing with a situation where everything is scaled on a
% range of -1 to 1
try is_wav(1); catch is_wav = 0; end

try nbits(1); catch nbits = 16; end

max_dB = 20*log10(2^nbits);
final_amp = 2^(nbits-1)/(10^((max_dB-dB)/20));

% Remove any dc-component
wvf = detrend(wvf,0);

% If we are dealing with a waveform where the maximum absolute value is <1 then
% scale to maximum representation before scaling to target dB value.
if all(max(abs(wvf))) <= 1
  IS_WAV = 1;
  wvf = wvf * 2^(nbits-1);
end

% Calculate the RMS of the input waveform
wvf_rms = sqrt(mean(wvf.^2));

% Calculate an RMS scaling factor
scale_factor = final_amp*sin(pi/4)./wvf_rms;

% Apply the scaling factor
scaled = wvf .* repmat(scale_factor,size(wvf,1),1);

if is_wav
  scaled = scaled ./ (2^(nbits-1));
end

return

