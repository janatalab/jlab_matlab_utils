function ev_info = create_ev_info
% ev_info = create_ev_info;
%
% Creates a basic explanatory variable(s) (EV) information structure

% 03/20/06 Petr Janata

ev_info.set_type = '';  % category name for this set of EVs
ev_info.ev_names = {''};
ev_info.nev = 0;
ev_info.fnames = {''};
