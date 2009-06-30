% template describing the 'ri' struct output from each physio reader

irun = 1;
channel_name = 'physio_channel';

% from read_keithley - all indices represent samples
rk.ri(irun)
rk.ri(irun).slice_onsets
rk.ri(irun).pos_events
rk.ri(irun).neg_events
rk.ri(irun).key_events
rk.ri(irun).cardiac
rk.ri(irun).respir

% from read_mate - vol onsets are in seconds, not samples
rm.ri(irun)
rm.ri(irun).vol_onsets
rm.ri(irun).cardiac
rm.ri(irun)

% from read_biopac - all indices are in samples, can be converted using
% .meta.srate or .meta.stime
rb.ri(irun)
rb.ri(irun).meta
rb.ri(irun).meta.stime
rb.ri(irun).meta.srate
rb.ri(irun).meta.baseline_samps
rb.ri(irun).signal
rb.ri(irun).signal.(channel_name)
rb.ri(irun).trigger.onsets
rb.ri(irun).trigger.offsets
