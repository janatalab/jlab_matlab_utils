function varargout = reson_filter(varargin)
%
% processes a signal through a reson filter
%
% INPUTS (passed in as tag,value pairs)
%   inputSig:      input signal to process with the filter
%   resonantFreq:  resonant frequency of the filter
%   bandwidth:     half-power bandwidth of the filter in Hz.
%   Fs:            sampling rate of the signal
%   gain:          constant (between 0 and 1) that is multiplied by
%                  numerator of transfer function. Default = 1.
%   resonType:     either 'reson' or 'reson_z'. reson_z adds two zeros
%                  to the filter in order to improve the low frequency rolloff
%                  (for low frequency filters), or high frequency rolloff (for
%                  high frequency filters). Default is 'reson_z' if omitted.
%   normalize:     logical 1 or 0 for true or false. Default = 1.
%   Qfactor:       Ratio of resonantFreq/bandwidth. Either Qfactor or
%                  bandwidth may be passed in as arguments. If both are passed in,
%                  make sure they are compatible (i.e. Q=Fr/bw)
%
%
% OUTPUTS
% If inputSig is passed in, outputs are [outputSig,B,A], where B and A are the filter coefficients of
% the numerator and denominator. If inputSig is not passed in, then
% the output simply contains the filter coefficients [B,A].
%
%
% 1/23/08 Stefan Tomic. Took code out of BTB into a more general
%                       purpose function
%
 
for iarg = 1:2:nargin
  switch(varargin{iarg})
   
   case 'inputSig'
    inSig = varargin{iarg+1};
  
   case 'resonantFreq'
    Fr = varargin{iarg+1};
    
   case 'bandwidth'
    bw = varargin{iarg+1};
    
   case 'gain'
    gain = varargin{iarg+1};
    
   case 'resonType'
    resonType = varargin{iarg+1};
    
   case 'normalize'
    normalize = varargin{iarg+1};
    
   case 'Fs'
    Fs = varargin{iarg+1};
    
   case 'Qfactor'
    Qfactor = varargin{iarg+1};
    
  end
end

try
  normalize;
catch
  normalize = 1;
end

try
  gain;
catch
  gain = 1;
end

if(resonantFreq > Fs/2)
  error('Resonant freqency is above Nyquist frequency');
end

if(gain < 0 || gain > 1)
  error('Gain is not between 0 and 1.');
end

%if Qfactor was passed instead of bw, calculate bw from Qfactor.
%if both Qfactor and bw were passed in, make sure that they are compatible.
%if only bw was passed in, then no need to calculate Qfactor, since
%we only need bw from here on.
if(~exist('bw','var') && exist('Qfactor','var'))
  bw = Fr/Qfactor;
elseif(exist('bw','var') && exist('Qfactor','var'))
  if(bw ~= Fr/Qfactor)
    error('Incompatible bandwidth and Qfactor were passed in');
  end
end
  
R = 1 - bw./Fs*pi;  

%psi is the desired resonant frequency. Note that with reson
%filters that this differs slightly from the location of the pole
%angle. So we need to calculate the pole angle theta from the
%desired resonant frequency angle, psi.
psi = Fr/Fs*(2*pi);    

if(normalize)
  %we will set a value for normFactor here, which will normalize
  %the filter's magnitude response. The gain (passed in as a param)
  %can also be imposed on the normalized response, so that we can
  %scale the magnitude response effectively. Since non-normalized
  %reson filters' peak magnitude response is dependent on the
  %resonant frequency, normalization gives us better control over
  %the magnitude response.
  
  switch(resonType)
    %normalization method depends on the type of reson filter used
   case 'reson_r'
    normFactor = 1-R; 
   case 'reson_z'
    normFactor = (1-R^2)/2;
   case 'reson'
    normFactor = (1-R^2)*sin(theta);
  end
else
  normFactor = 1;
end

switch(resonType)
 case 'reson_z'
  B_unscaled = [1 0 -1];
  % approximating actual pole angle theta from desired peak response
  % (Steiglitz, 1994)

  theta = acos( (1+R^2)./(2*R).*cos(psi) );
 case 'reson'
  B_unscaled = [1];
  theta = acos( (2*R)./(1+R^2).*cos(psi) );
 case 'reson_r'
  psi = theta; %we have no peak response correction equation for reson_r
  B_unscaled = [1 0 -R];
end


B = B_unscaled.*gain.*normFactor;
A = [1 -2*R*cos(theta) R^2];

if(exist('inSig','var'))
  varargout{1} = filter(B,A,inSig);   
  varargout{2} = B;
  varargout{3} = A;
else
  varargout{1} = B;
  varargout{2} = A;
end
