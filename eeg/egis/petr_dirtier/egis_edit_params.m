function ep = edit_params
%
%

ep.start_samp = 1;
ep.stop_samp = [];
ep.do_eyemove = 0;

%
% Define defaults for adult 129 channel GSN
%

ep.inflfteye = 127;
ep.infrteye = 126;
ep.suplfteye = 26;
ep.suprteye = 8;
ep.horizlfteye = 128;
ep.horizrteye = 125;

%
% Work on all cells by default
%
ep.do_cells = [];		% Empty does all cells

%
% Default threshold parameters.  These are defined in terms of standard
% deviation in the distributions that are computed from the data.
%

ep.thresh.blink = 2.5;
ep.thresh.transient = 2.5;
ep.thresh.absvolt = 2.5;

ep.thresh.numstd = 2.5;  % exclude channel if any value
                               % exceeds the number of standard deviations 
                               % specified.  0 means no exclusion

ep.thresh.stdmag = 2.5;  % Reject channel if std dev on it exceeds this
                                % value

%ep = class(ep,'egis_edit_params');
