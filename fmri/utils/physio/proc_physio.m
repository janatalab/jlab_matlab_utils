function phys = proc_physio(ps)
% phys = proc_physio(subid, expanme);
%
% ps - is a structure returned by init_phys_struct which contains all of the
%      information necessary to perform the conversion.
%
% Adapted from modfmri/convert_physio.m
%
% Wrapper for script which converts ASCII file with physiological monitoring
% data into a mat file that can then be used for design matrix creation.
%
% The intent of this function is for it to ultimately handle multiple types of
% input data formats
%

% 10/27/05 Petr Janata - adapted script
insane = 0;  % Flag for something having gone quite wrong
phys = [];

% Extract some variables from the ps struct
sinfo = ps.sinfo;
sid = sinfo.id;
pp = ps.pp;
logfid = ps.logfid;

phys_fname = fullfile(ps.datapath, sid, sprintf('%s%s%s', sid,ps.phys_data_stub,ps.phys_data_suffix));
check_exist(phys_fname);

% Call the parser
ri = read_physmon(phys_fname, ps.pp);

% Eliminate a run up front if desired
if ~isempty(ps.elim_runs)
  fprintf(logfid, 'Eliminating run(s) %s\n', sprintf('%d ', ps.elim_runs));
  ri(ps.elim_runs) = [];
end

% Number of runs is determined by runs flagged for analysis in
% get_<expname>_sinfo.m
use_runs = sinfo.use_runs;
nruns = length(use_runs);

if nruns ~= length(ri)
  error_str = sprintf('ERROR:proc_physio:Expected %d runs, detected %d\n', nruns, length(ri));
  fprintf(logfid, error_str);
  error(error_str);
  return
else
  for irun = 1:nruns
    ri(irun).id = use_runs(irun);
  end
end

%
% NOTE, from here on out, we assume that the entries in ri reflect the runs
% specified in sinfo.use_run 
%

%
% Perform some data integrity checks
%
fprintf(logfid, 'Performing data integrity checks ...\n');


ps.expect.nslice = sum(sinfo.nvol(use_runs)*ps.nslice_per_vol);

if nruns ~= ps.expect.nruns
  fprintf(logfid,'WARNING:proc_physio:Really wanted %d runs; found %d\n', ps.expect.nruns, nruns);
  insane = 0;
  
  % Change expectations for number of events given non-canonical number of
  % runs.
  ps.expect.pos_events = nruns*ps.expect.pos_events_per_run;
  ps.expect.neg_events = nruns*ps.expect.neg_events_per_run;
end

npos = length(cat(1,ri.pos_events));
if npos ~= ps.expect.pos_events
  fprintf(logfid,'WARNING:proc_physio:Expected %d positive markers, detected %d\n', ps.expect.pos_events, npos);
  insane = 1;
end

nneg = length(cat(1,ri.neg_events));
if ~any(ismember(nneg, ps.expect.neg_events))
  fprintf(logfid,'WARNING:proc_physio:Expected %d negative markers, detected %d\n', ps.expect.neg_events(1), nneg);
  insane = 1;
end

nslice = length(cat(1,ri.slice_onsets));
nslice_per_vol = ps.nslice_per_vol;
if nslice ~= ps.expect.nslice
  for irun = 1:nruns
    run_idx = ri(irun).id;  % use this to index into original sinfo structures
    expect_nslice = sinfo.nvol(run_idx)*nslice_per_vol;
    actual_nslice = length(ri(irun).slice_onsets);
    if actual_nslice ~= expect_nslice
      fprintf(logfid,'WARNING:proc_physio:Run %d: Expected %d slices, detected %d\n', irun, expect_nslice, actual_nslice);

      if expect_nslice == actual_nslice-nslice_per_vol
	ndiff = actual_nslice-expect_nslice;
	fprintf(logfid,'Eliminating %d excess slices from beginning of run %d\n', ndiff, irun);
	
	% When volumes are dropped, they are dropped at the beginning of the
	% run.  Therefore, all timing info has to be changed to accommodate this.
	ri(irun).slice_onsets(1:ndiff) = [];
	
	new_start_time = ri(irun).slice_onsets(1);
	
	ri(irun).slice_onsets = ri(irun).slice_onsets - new_start_time;
	
	ri(irun).pos_events(2:end) = ...
	    ri(irun).pos_events(2:end)-new_start_time;
	ri(irun).neg_events = ri(irun).neg_events - new_start_time;
	ri(irun).key_events = ri(irun).key_events - new_start_time;
	
	ri(irun).cardiac = ri(irun).cardiac - new_start_time;
	elim_idx = find(ri(irun).cardiac < 0);
	ri(irun).cardiac(elim_idx) = [];
	
	elim_idx = find(((0:length(ri(irun).respir)-1)/pp.Fs ...
	    -new_start_time) < 0);
	ri(irun).respir(elim_idx) = [];
      elseif expect_nslice < actual_nslice
	% Scanner crash scenario
	ndiff = actual_nslice-expect_nslice;
	fprintf(logfid,'Eliminating %d excess slices from end of run %d\n', ndiff, irun);
	ri(irun).slice_onsets(actual_nslice ...
	    -(actual_nslice-expect_nslice)+1:end) = []; 
	
	% Eliminate all responses and event markers that occurred after
	% scanner quit saving the data
	stop_time = ri(irun).slice_onsets(end);
	bad_idx = find(ri(irun).pos_events > stop_time);
	ri(irun).pos_events(bad_idx) = [];

	bad_idx = find(ri(irun).neg_events > stop_time);
	ri(irun).neg_events(bad_idx) = [];

	bad_idx = find(ri(irun).key_events > stop_time);
	ri(irun).key_events(bad_idx) = [];

	% Let us capture next heartbeat so that we can compute phase in
	% cardiac cycle during design matrix construction
	bad_idx = find(ri(irun).cardiac > stop_time + ps.cardiac_slop_s);
	ri(irun).cardiac(bad_idx) = [];

	bad_idx = find((0:length(ri(irun).respir)-1)/pp.Fs > stop_time);
	ri(irun).respir(bad_idx) = [];
	
      else
	insane = 1;
      end
    end
  end % for irun = 1:nruns
end % if nslice ~= expect.nslice

ps.insane = insane;
if insane
  fprintf(logfid,'\tErrors detected!\n');
else
  fprintf(logfid,'\tNo errors found\n');
end

% Pack the return structure
phys.ps = ps;
phys.ri = ri;

if pp.savedata
  phys_matname = fullfile(ps.datapath, sid, sprintf('%s_phys.mat', sid));
  fprintf(logfid,'Saving physio information to: %s\n', phys_matname);
  save(phys_matname, 'phys')
end

return