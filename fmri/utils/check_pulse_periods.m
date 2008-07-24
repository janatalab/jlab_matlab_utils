function [r] = check_pulse_periods(pulse_times,params,r)
% [r] = check_pulse_periods(pulse_times,params)
%
% Given a vector of event onset times, e.g. MRI scanner pulses, that are
% expected to have a constant period, the script identifies aberrent pulse periods
%
% pulse_times - vector of pulse times (msec)
%
% r is the reporting structure which contains a data structure with the
% following fields:
% bad_pulse_period_mask - logical vector indicating which pulse periods are
%                         aberrant.  The first element in this vector refers to
%                         the period between the first and second elements in
%                         the pulse_times vector.
% pulse_period_too_short_mask - those periods that were too short
% pulse_period_too_long_mask - those periods that were too long
% nmiss - estimate of number of missed pulses
% median_period - estimate of the median pulse period
% msglog - list of all messages that were generated

% 05/12/06 Petr Janata - updated code with additional error checking
% 08/03/06 PJ - implemented within the emerging reporting framework 
% 12/11/07 PJ - fixed minor bug in which nmiss was not being created when
%               no pulses were missed

r.type = 'pulse_period_check';  % Identify the type of this reporting instance
r.report_on_fly = 1;

try max_pulse_dev = params.max_pulse_dev; catch
  max_pulse_dev = 20; % maximum deviation (in msec) of a pulse from the expected
                     % period.
end

try min_pulse_period = params.min_pulse_period; catch
  min_pulse_period = 50;  % won't be going faster than 20 slices/sec
end

% Look at variability of timing pulses, and check for skipped pulses
pulse_periods = diff(pulse_times);  % determine all of the pulse periods
median_period = median(pulse_periods);

% Check to see if the median pulse period is too short
pulse_period_too_short_mask = zeros(size(pulse_periods));
if median_period < min_pulse_period
  msg = sprintf(['Median pulse period (%4.1f) is less than minimum allowed ' ...
	'period (%4.1f)'], median_period, min_pulse_period)
  r = update_report(r,msg);
 
  pulse_period_too_short_mask = pulse_periods < min_pulse_period;
  msg = sprintf('Detected %d pulse periods that were too short\n', sum(pulse_period_too_short_mask));
  r = update_report(r,msg);
  
  % Revise the median pulse period estimate
  median_period = median(pulse_periods(~pulse_period_too_short_mask));
end

% Check for pulse periods that are too long
pulse_period_too_long_mask = pulse_periods > median_period+max_pulse_dev;
bad_pulse_period_idxs = find(pulse_period_too_long_mask);

nbad = length(bad_pulse_period_idxs);
nmiss = 0;
if nbad
  msg = sprintf('Detected %d pulse periods that were too long\n', nbad);
  r = update_report(r,msg);

  % Estimate the number of missed pulses
  for ibad = 1:nbad
    nmiss = nmiss + round(pulse_periods(bad_pulse_period_idxs(ibad))/median_period)-1;
  end % for ibad=
  msg = sprintf('Estimate of %d missed pulses\n', nmiss);
  r = update_report(r,msg);
end

% Pulse periods that are too long aren't technically bad in that there are only
% pulses missing.  Short pulse periods are bad though, so add them to the bad_pulse_period_mask.
r.data.bad_pulse_period_mask = pulse_period_too_short_mask; % | pulse_period_too_long_mask;
r.data.pulse_period_too_short_mask = pulse_period_too_short_mask;
r.data.pulse_period_too_long_mask = pulse_period_too_long_mask;
r.data.nmiss = nmiss;
r.data.median_period = median_period;
