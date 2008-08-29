function new_trialdata = replace_trial(trialdata,good_chan,xelec,yelec,zelec)
% new_trialdata = replace_trial(trialdata,good_chan,xelec,yelec,zelec)
%
%eplace bad channels in a trial of data 
%		
%
% trialdata = data with channels as columns
%good_chan = good channel list
%xelec,yelec,zelec = x,y,z, electrode positions
if nargin ~= 5
	error('incorrect number of input arguments')
end
v = zeros(1,size(good_chan,2));
x = zeros(1,size(good_chan,2));
y = zeros(1,size(good_chan,2));
z = zeros(1,size(good_chan,2));
x = xelec(good_chan);
y = yelec(good_chan);
z = zelec(good_chan);
welec= 1;
[k,kinv,a,ainv,e]= k_and_e(welec,x,y,z);
for samp = 1:size(trialdata,1)
	w = trialdata(samp,:);
	v = w(good_chan);
	[p,q,error_check]= mateqs(welec,x,y,z,v,k,kinv,a,ainv,e);
	v_lap = interp_3d(welec,x,y,z,xelec,yelec,zelec,p,q);
	new_trialdata(samp,:) = v_lap;
end; %for samp = 1:chdr(c,NPoints)

