function p = init_physmon_params(config)
% params = init_physmon_params(config);
%
% Initializes parameter structure for reading physiological monitoring files.
%
% config - specifies a particular physiological monitoring system/configuration.
%
% See also README_physio.txt

% 1/10/02 PJ
% 05/27/04 PJ  Added .is_event_trig_pos and .is_button_trig_pos fields
% 10/25/04 PJ  Added .neg_trig_thresh_key
% 07/20/06 PJ  Modified to handle multiple acquisition systems

p.loaddata = 1;
p.savedata = 0;
p.outsuffix = '';
p.savefixed = 0;			% save data after it has been adjusted
                                        % for dropped volumes, etc.

switch config
  case 'keithley'
    p.Fs = 250;				% sampling rate (Hz)
    p.bipolar = 0;			% are event marker pulses bipolar?
    p.pos_trig_thresh = 2;
    p.neg_trig_thresh = -2;
    p.neg_trig_thresh_key = -4;		% threshold for key press events
    p.is_event_trig_pos = 1;
    p.is_button_trig_pos = 1;
    
    % if markers are bipolar, this is the window within which negative/positive transition should occur    
    p.event_window_ms = 30;			
    
    % time from EPI start to first stimulus marker
    p.epi_start_to_first_stim = 17500;	
    
    % tolerance for marker locations
    p.event_slop = 2;		
    
    % for finding runs based on receiver unblank maximum duration of a run, in seconds)
    p.slice_diff_thresh = 100;		

    % Channel assignments
    p.cardiac_chan = 1;
    p.respir_chan = 2;
    p.slice_chan = 3;
    p.button_chan = 4;
    p.event_chan = 5;

    % Removal of mean from event channels prior to locating trigger evenets
    p.remove_dc = 1;		

    % Slice timing channel parameters
    p.TR = 3;
    p.nslice_per_vol = 27;			% number of slices per volume

    % number of dummy volumes at the beginning of each run
    p.ndummy_vol = 2;			

    % Physiological channel post-processing

    % cutoff frequency for filtering respiratory data. Leave empty if no filtering is desired
    p.low_pass_cutoff = 5;		
					
  case 'mate'
    p.source = 'mate';
    
    % Scanner trigger related variables
    p.trig_thresh = -2.5;
    p.trig_dir = 'neg';
    p.Fs_trig = 100;
    p.trig_slop_sec = 0.02; 
    p.trig_is_vol = 1;  % volume or slice triggers
    
    p.pulse_thresh = -2.5;
    p.pulse_dir = 'neg';
    p.Fs_pulse = 50;
    p.Fs_resp = 50;

  case 'biopac'
    p.trig_thres = 2.5;
    p.trig_dir = 'pos';
    p.trig_is_vol = 1;

    % Channel assignments
    p.scr_chan = 1;
    p.respir_chan = 2;
    p.cardiac_chan = 4;
    p.slice_chan = 7;
    p.block_chan = 9;
    p.event_chan = 14;
    
  otherwise
    fprintf('Unknown configuration: %s\n', config)
end
