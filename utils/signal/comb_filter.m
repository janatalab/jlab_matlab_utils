function varargout = comb_filter(varargin)
%
% processes a signal through a comb filter
%
% INPUTS (passed in as tag,value pairs)
%   inputSig:       input signal to process with the filter
%   resonantFreq:   lowest resonant frequency of the filter
%                   (fundamental frequency of the comb filter)
%   Fs:             sampling rate of the signal
%   halfEnergyTime: Only configured for feedback comb filter (combType='normal').
%                   Tunes the filter so that it has a  half energy
%                   decay time specified in seconds (see Scheirer, 1998)
%   combType:       either 'normal' or 'inverse'. 'inverse'
%                   implements a feedfoward comb filter with a magnitude response
%                   that is the inverse of the feedback ('normal') version. Default
%                   is 'normal'.
%   gain:           constant (between 0 and 1) that is multiplied by
%                   numerator of transfer function. Default = 1.
%   poleRadius:     If 'halfEnergyTime' is set, this argument is
%                   ignored. Radius of the poles for the feedback comb filter (combType='normal').
%                   This affects the overall gain of the comb filter.
%   zeroRadius:     Radius of the zeros for the feedforward comb filter (combType='inverse')
%                   This affects the overall gain of the comb filter.
%
% OUTPUTS
% If inputSig is passed in, outputs are [outputSig,B,A], where B and A are the filter coefficients of
% the numerator and denominator. If inputSig is not passed in, then
% the output simply contains the filter coefficients [B,A].
%
%
% 1/26/08 Stefan Tomic. Took code out of BTB into a more general
%                       purpose function
 
for iarg = 1:2:nargin
  switch(varargin{iarg})
   
   case 'inputSig'
    inSig = varargin{iarg+1};
  
   case 'resonantFreq'
    Fr = varargin{iarg+1};
    
   case 'gain'
    gain = varargin{iarg+1};
    
   case 'halfEnergyTime'
    tHalf = varargin{iarg+1};
    
   case 'combType'
    combType = varargin{iarg+1};
    
   case 'Fs'
    Fs = varargin{iarg+1};
    
   case 'poleRadius'
    poleRadius = varargin{iarg+1};
    
   case 'zeroRadius'
    zeroRadius = varargin{iarg+1};
    
  end
end

try
  combType;
catch
  combType = 'normal';
end

try
  gain;
catch
  gain = 1;
end
     
switch(combType)
  
  case 'normal'

   if(exist('tHalf','var') && exist('poleRadius','var'))
     warning(['Both half energy decay time and pole radius were specified.' ...
	      ' Ignoring the pole radius value.']);
   end

   Tr = 1./Fr;
   trSamps = Tr*Fs;

   if(exist('tHalf','var'))
     tHalfSamps = tHalf*Fs;
     alpha = 0.5.^(trSamps./tHalfSamps);
     A = [1 zeros(1,trSamps-1) -1*alpha];
     B = 1-alpha;
   elseif(exist('poleRadius','var'))
     A = [1 zeros(1,trSamps-1) -1*poleRadius^trSamps];
     B = gain;
   else
     error('comb_filter needs either a half energy time or pole radius');
   end

 case 'inverse'
  
  %for the inverse comb filter, the first zero is actually at twice
  %the angle of the lowest resonance, so the delay is half of the
  %period corresponding to the first resonance.
  Tr = 1./Fr;
  trSamps = Tr*Fs./2;
  
  B =  gain.*[1 zeros(1,trSamps-1) -1*zeroRadius^trSamps];
  A = 1;
  
end
    
varargout = {};
if(exist('inSig','var'))
  varargout{end+1} = filter(B,A,inSig);
end

varargout{end+1} = B;
varargout{end+1} = A;
    
return
