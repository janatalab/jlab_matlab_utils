function avref_trialdata= average_reference(trialdata,mask)
%avref_trialdata= average_reference(trialdata,mask)
%
%average references trial
%

%Modification history
%
% written 7/17/95 by RS
%
% cleaned out a couple of loops and messed with it in the process -- PJ 7/18/99
% resized to accommodate new input size of mask - RS 
%
% 7/15/01 PJ Cleaned up
[npoints NChan] = size(trialdata);

trialdata(:,NChan+1) = zeros(npoints,1);

temptrialdata = zeros(npoints,NChan+1);
temptrialdata =trialdata.* repmat(mask,npoints,1);
divisor=sum(mask);
average=sum(temptrialdata')/divisor;
channel=-average';

avref_trialdata = temptrialdata + repmat(channel,1,NChan+1); 
avref_trialdata=avref_trialdata .* repmat(mask,npoints,1);


