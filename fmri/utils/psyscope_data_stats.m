function psyscope_data_stats(data, verbosity)
%
% psyscope_data_stats(data, verbosity);
%
%  Assumes data structure returned by rd_psyscope_data()
%

% 03/24/00 PJ

nevents = length(data.trial);

disp(sprintf('%d events', nevents))

trials = unique(data.trial);
disp(sprintf('Encountered %d unique trials', length(trials)))

states = unique(data.state);
disp(sprintf('Encountered %d unique state_masks', length(states)))

if verbosity > 1
  disp(sprintf('\t%d', states))
end

%
%  Determine the number of active input lines
%

active =find(states < 512);

mask = 0;

for i = 1:length(active)
  mask = bitor(mask, states(active(i)));
end

disp(sprintf('Active lines: %s', sprintf('%d ', find(bitget(mask,1:8)))))
