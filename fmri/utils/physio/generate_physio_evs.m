function phys = generate_physio_evs(phys)
% phys = generate_physio_evs(phys);
%
% Generates explanatory variables (EVs/regressors) with physiological data.
% These regressors can then be loaded into regression models for fMRI data.
%
% phys - a structure containing all of the relevant data and metadata
%
% See also: populate_phys_struct, proc_physio

% 11/01/05 Petr Janata

% copy some often used values from the structure
logfid = phys.ps.logfid;

subid = phys.ps.sinfo.id;
ri = phys.ri;
pp = phys.ps.pp;

nruns = length(ri);
nslice_per_vol = phys.ps.pp.nslice_per_vol;

% Determine how many sets of regressors we are generating
nsets = length(phys.out_types);

% Loop over all of the runs
for irun = 1:nruns
  % We need to refer to runs in terms of their original indices, i.e. order in
  % which they were acquired.
  runid = phys.ps.sinfo.use_runs(irun); 
  
  nslices = length(phys.ri(irun).slice_onsets);
  nvol = phys.ps.sinfo.nvol(runid);

  % Perform a sanity check
  if nslices ~= nvol*nslice_per_vol
    fprintf(logfid,'ERROR:generate_physio_evs:Number of slices (%d) does not match expected number of slices (%d)', nslices, nvol*nslice_per_vol);
    error
  end
               
  ev_info = [];
  
  % Loop over all of the desired regressor sets
  for iset = 1:nsets
    set_type = phys.out_types{iset};
    outstub = sprintf('%s_run%d_%s', subid, runid, set_type);

    % Keep track of EV info in the ri (run info) structure
    ev_info(iset).set_type = set_type;

    % Do any necessary preprocessing and setting up of variables
    switch set_type
      case 'cardiac'
	curr_regress = zeros(nslices,6);
	ev_info(iset).ev_names = {'Sin1','Sin2','Sin3','Cos1','Cos2','Cos3'};
      case {'respir', 'respiration'}
	curr_regress = zeros(nslices,1);
	ev_info(iset).ev_names = {'Respiration'};
      otherwise
    end
    
    % Loop over all of the slices in the volume
    for islice = 1:nslices
      % Determine the sample offset of this slice on the original physiological
      % data timescale
      samp_idx = fix(ri(irun).slice_onsets(islice)*pp.Fs)+1;

      switch set_type
	case 'cardiac'
	  card1_idx = max(find(ri(irun).cardiac <= ri(irun).slice_onsets(islice)));
	  if isempty(card1_idx)
	    card1 = 0;
	  else
	    card1 = ri(irun).cardiac(card1_idx);
	  end
	  card2_idx = min(find(ri(irun).cardiac > ...
	      ri(irun).slice_onsets(islice)));

	  % If physio collecting is turned off part way through last cardiac
	  % cycle, we need to estimate when in the last cardiac cycle the slice occurred.
	  if isempty(card2_idx)
	    next_card = ri(irun).cardiac(end) + ...
		mean(diff(ri(irun).cardiac(end-9:end)));
	    % Make sure nothing is really wrong
	    if ri(irun).slice_onsets(end) > next_card
	      error('Something is wrong with cardiac data')
	    else
	      card2 = next_card;
	    end
	  else
	    card2 = ri(irun).cardiac(card2_idx);
	  end
	  % As described in Dagli et al. 1999, NeuroImage 9, 407-415
	  theta = (ri(irun).slice_onsets(islice) - card1)/(card2-card1);
	  curr_regress(islice,:) = [sin(2*pi*theta*(1:3)) cos(2*pi*theta*(1:3))];
	  
	case {'respir', 'respiration'}
	  % Figure out the range of samples in the respiration waveform (at the
	  % sampling rate of the physio file) that we need to consider with
	  % regard to the acquisition of this particular slice
	  start_samp = samp_idx-phys.ps.nsamp_resp_peri;
	  if start_samp < 1
	    start_samp = 1;
	  end
	  stop_samp = samp_idx+phys.ps.nsamp_resp_peri;
	  if stop_samp > length(ri(irun).respir)
	    stop_samp = length(ri(irun).respir);
	  end
	  
	  % Calculate the value for this slice
	  curr_regress(islice) = mean(ri(irun).respir(start_samp:stop_samp));
	  
	otherwise
	  fprintf(logfid, 'Do not know how to handle regressor set: %s\n', set_type);
      end % switch phys.out_types{iset}
    end % for islice=

    % Post-process the regressors
    switch set_type
      case {'respir', 'respiration'}
	curr_regress = detrend(curr_regress);
      otherwise
    end
  
    % Rearrange into a 3D design matrix which is of the form 
    % nvol X nev X nslice
    regset = reshape(curr_regress, nslice_per_vol, nvol, size(curr_regress,2));
    
    regset = permute(regset,[2 3 1]);
    
    % Correct for slice acquisition order
    regset(:,:,phys.ps.slice_order) = regset;
    nev = size(regset,2);
    ev_info(iset).nev = nev;
  
    for islice = 1:nslice_per_vol
      for iev = 1:nev
	outfname = fullfile(phys.outdir, sprintf('%s_ev%d_slice%02d.txt', outstub, iev, islice));
	fprintf(logfid,'Saving EV to file: %s\n', outfname);
	
	ev_info(iset).fnames{islice,iev} = outfname;
	
	outdata = regset(:,iev,islice);
	switch phys.file_format
	  case 'fsl'
	    save(outfname,'outdata','-ascii')
	  otherwise
	    fprintf(logfid,'Do not know how to handle desired EV output format: %s\n', phys.file_format);
	end % switch phys.file_format
      end % for iev
    end % for islice
  end % for iset=
  
  phys.ri(irun).ev_info = ev_info;
end % for irun=