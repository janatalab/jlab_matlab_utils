function vol_trialdata = cal_gain(trialdata,CalGain,CalZero)
%vol_trialdata = cal_gain(trialdata,CalGain,CalZero)
%
%corrects for calibration zeros and gains returning microvolts

%MODIFICATION HISTORY
%
%	written 7/17/99 R.S.
%
%	modified 7/19/95 R.S. and (P.J.)
%   removing bad FORTRAN related coding habits.

size_trial = size(trialdata);
for i=1:size_trial(2)
	CalGain(i) = 38/0.4883;
	end;
	
CalZeros_mat = ones(size_trial(1),1)*CalZero';
CalGain_mat = ones(size_trial(1),1)*CalGain';

vol_trialdata = (trialdata - CalZeros_mat).*38./CalGain_mat;



