function ps = init_phys_struct(config)
% ps = init_phys_struct()
%
% Initializes a structure with various control variables for reading files
% containing accessory physiological data such as heart rate, respiration, etc.
%
% Initially written to support analyses of fMRI data
%

% 10/28/05 Petr Janata
% 03/16/08 FB - added config, passed on to init_physmon_params

ps.dataroot = '';  % common global root of all data files
ps.datapath = '';  % experiment specific path
ps.phys_data_stub = '';
ps.phys_data_suffix = '';
ps.logfid = 1;

ps.cardiac_slop_s = 2;
ps.nslice_per_vol = [];
ps.elim_runs = [];

% Initialize more specific parsing parameters
ps.pp = init_physmon_params(config);

expect.nruns = [];
expect.pos_events_per_run = [];
expect.pos_events = [];
expect.neg_events_per_run = [];
expect.neg_events = [];

ps.expect = expect;

ps.verbose = 0;