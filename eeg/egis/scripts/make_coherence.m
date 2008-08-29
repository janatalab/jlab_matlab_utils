function [coh_mat] = make_coherence(csdmfiles,outfname,win, cf_list, cf_flag, bad_chans);
%[coh_mat] = make_coherence(csdmfiles,outfname,win, cf_list, cf_flag, bad_chans);
%
%  this function makes coherence files derived files from csdm files
%
%  csdmfiles = the names of the csdm files must be stored in an array where each row 
%              is a filename.   
%  outfname = root name of output file.  The number of the channel or frequency with
%  respect to which the coherence matrix is computed is appended to the
%  root name.  _f3 indicates frequency bin 3 of the original csdm.
%	_c69 indicates channel 69.
%			
% win = 	which window (observation number) to use in each file 
%		if win == [], all observations are used
%		if length(win) > 0 and there is more than one csdmfile specified, a
%           separate file is created for each window, with NObs in each file
%			corresponding to the number of input csdmfiles
%
%	cf_list = list of channels or frequency bins
%	cf_flag = 'chan' if coherence is computed relative to one channel.  Resulting matrix
% 				is Npoints X NChan
%		 	= 'freq' if coherence is computed relative to all channels, but only for a 
%				single frequency.  Resulting matrix is NChan X NChan
%
%	bad_chans = list of bad channels -- one row for each row of csdmfiles

%  01/20/96 PJ Started modification tracking and modified for propagation of bad_chans
%  01/27/96 PJ Combined make_coh_by_freq and make_coh_by_chan

status = -1;
format compact

%
%  Check input arguments
%

if (nargin<5)
  error('Too few arguments to make_coherence()')
elseif nargin<6
  bad_chans = [];
else
  if size(bad_chans,1) ~= size(csdmfiles,1)
    error('csdmfile list and bad_chans lists don''t match');
  end
end;

if ~(strcmp(cf_flag, 'freq') | strcmp(cf_flag, 'chan'))
  error(['Invalid cf_flag: ' cf_flag]);
end

ave_hdr_offsets_v;

%
% Read in copy of the csdm header
%

[csdm_fhdr,csdm_chdr]=get_csdm_hdr(csdmfiles(1,:));

if win == []
  win = 1:csdm_chdr(1,NObs);
end


%
% Copy csdm headers to average header
%

[ave_fhdr,ave_chdr,ave_ename,ave_cnames,ave_fcom,ave_ftext] = csdmhdr2avehdr(csdmfiles, win);

% Correct the NPoints field for each cell header
% This must be done independently in each m-file since in some cases the 
% resulting matrices are NChan X NChan instead of NPoints X NChan

for c = 1:csdm_fhdr(NCells)
  if strcmp(cf_flag, 'chan')
	ave_chdr(c,NPoints) = csdm_chdr(c,NPoints)/2;
  else
 	ave_chdr(c,NPoints) = ave_fhdr(NChan);
  end 
end

%
% Determine the number of output files corresponding to windows
%

if size(csdmfiles,1) > 1 & length(win) > 1
  nwin_files = length(win);
else
  nwin_files = 1;
end

nitems = size(cf_list,2);
nwin = length(win);

%
% Get fileids for all of the resulting coherence files
% The file ids are stored in a nitems X nwin matrix
%

new_root = outfname;
if strcmp(cf_flag, 'freq')
  new_root = [new_root '_cf'];  % cf = coherence by frequency
else
  new_root = [new_root '_cc'];	% cc = coherence by channel
end

fid_mat = zeros(nitems,nwin_files);

for iref = 1:nitems
	iref_str = int2str(cf_list(iref));
	outfname1 = [new_root iref_str];
	if nwin_files == 1
	  fid_mat(iref, 1) = fopen(outfname1,'wb');
	  wt_ave_hdr_v(fid_mat(iref,1),ave_fhdr,ave_chdr,ave_ename,ave_cnames,ave_fcom,ave_ftext);
	  fseek(fid_mat(iref, 1), 0, 'eof');
	else
	  for fref = 1:nwin_files
	    text_ref = int2str(fref);
	    outfname1 = [new_root iref_str '_w' text_ref];
	    fid_mat(iref, fref) = fopen(outfname1,'wb');
	    wt_ave_hdr_v(fid_mat(iref,fref),ave_fhdr,ave_chdr,ave_ename,ave_cnames,ave_fcom,ave_ftext);
%	    fseek(fid_mat(iref, fref), 0, 'eof');
		fclose(fid_mat(iref,fref));
	  end
	end
end;
ch_pair_indices;

ch_str = [];
for c = 1:length(cf_list)
  ch_str = [ch_str int2str(cf_list(c)) ' '];
end
if strcmp(cf_flag, 'freq')
  disp(['Computing all channel coherences for frequency bins ' ch_str]);
else
  disp(['Computing coherences for all frequency bins with respect to chans: ' ch_str]);
end

[cell_data_offsets]=get_csdm_cell_off(csdm_fhdr, csdm_chdr);

for icell = 1:csdm_fhdr(NCells)
  disp(['Working on cell ' int2str(icell)]);
  npts = ave_chdr(icell,NPoints);
  nchan = ave_fhdr(NChan);
	
  for ifile = 1:size(csdmfiles,1)
    disp(['  Input file: ' csdmfiles(ifile,:)]);
    for iobs = 1:nwin
      disp(['  Observation ' int2str(iobs) '  (' int2str(win(iobs)) ')']);

		% read in the csdm matrix
	    [csdm_data, csdm_fhdr, csdm_chdr] = read_csdm(csdmfiles(ifile,:), icell, win(iobs));
		
		% compute coherence matrix successively for each item in cf_list
		for iref = 1:nitems
		  if strcmp(cf_flag, 'freq')
		    disp(['    Frequency band ' int2str(cf_list(iref))]);
		  else
		    disp(['    Channel ' int2str(cf_list(iref))]);
		  end
		  
		  % Size matrices
		  if strcmp(cf_flag, 'freq')
		    power_data = zeros(1,nchan);
		  else
		    power_data = zeros(npts, nchan);
		  end
		  
		  coh_data = zeros(npts,nchan);
		  numerator = zeros(npts, nchan);
	
		  % Generate the matrix of autospectra	
		  if strcmp(cf_flag, 'freq')					
		  	for ichan = 1:nchan
		  	  power_data(ichan) = csdm_data(cf_list(iref),chpair(ichan,ichan));
		    end;		  

		 	% Put in NaNs where there are bad_chans -- will be zeroed later
			if bad_chans ~= []
		 		power_data(bad_chans(ifile,find(bad_chans(ifile,:)))) = ones(size(find(bad_chans(ifile,:)))) * NaN;
		  	end
		  	power_data = power_data' * power_data;
		  else
		    for ichan = 1:nchan 
			  power_data(:,ichan) = csdm_data(1:npts,chpair(ichan,ichan));
			end
			if bad_chans ~= []
				power_data(:,bad_chans(ifile,find(bad_chans(ifile,:)))) = ones(size(power_data(:,bad_chans(ifile,find(bad_chans(ifile,:)))))) * NaN;
			end;
		  end
		  
		  % Make sure the autospectra matrix is valid
		  testreal = isreal(sum(power_data));
		  if testreal ~= 1
	        power_data(1,[6 54 96])
			error('complex numbers in power data matrix')
		  end;
			
		  % 
		  % PERFORM THE COHERENCE CALCULATION
		  %
		  
		  % COHERENCE BY FREQUENCY CASE
		  if strcmp(cf_flag, 'freq')
		 	% Generate the numerator for the coherence calc
		  	if strcmp(computer, 'MAC2')
		   	 for rchan = 1:nchan
			  numerator(rchan,:) = ppc_cmplx_mult(csdm_data(cf_list(iref),chpair(rchan, 1:nchan)),conj(csdm_data(cf_list(iref),chpair(rchan, 1:nchan))));
	         end
		    else
		      for rchan = 1:nchan
			numerator(rchan,:) = csdm_data(cf_list(iref),chpair(rchan, 1:nchan)).*conj(csdm_data(cf_list(iref),chpair(rchan, 1:nchan)));
				end
		  	end
		  
		  	% Compute the coherence
		 	coh_data = numerator./power_data;
		
		  else	% COHERENCE BY CHANNEL CASE
		  	% Generate the numerator for the coherence calc
			if strcmp(computer, 'MAC2')
			  for ichan = 1:nchan
			    numerator(:,ichan) = ppc_cmplx_mult(csdm_data(1:npts,chpair(cf_list(iref), ichan)),conj(csdm_data(1:npts,chpair(cf_list(iref), ichan))));
			  end
			else
			  for ichan = 1:nchan
			    numerator(:,ichan) = csdm_data(1:npts,chpair(cf_list(iref), ichan)).*conj(csdm_data(1:npts,chpair(cf_list(iref), ichan)));
			  end
			end
			
			% Compute the coherence
			coh_data = numerator./((power_data(1:npts,cf_list(iref))*ones(1,nchan)) .* power_data);
		  end
		  
		  %
		  % END OF COHERENCE CALC
		  %
		  
		  % Replace NaNs with zeros
		  coh_data(find(isnan(coh_data))) = zeros(size(find(isnan(coh_data))));
		  
		  %
		  % Error check various aspects of the coherence matrix
		  %
		  if strcmp(cf_flag, 'freq')
		  	testcoh = sum(diag(coh_data));
			if testcoh ~= (nchan - length(find(bad_chans(ifile,:))))
				if ~isempty(findstr(csdmfiles(ifile,:), '2D'))
				disp(['You are a moron who failed to detect bad chans' int2str(testcoh)])
			  	else
				disp(['Undetected Bad Chan ' int2str(testcoh)]);
				end
			end
		  else
			testcoh = sum(coh_data(1:npts,cf_list(iref)));
			if testcoh ~= npts
			  disp(['Your character is suspect ' int2str(testcoh)]);
			end
		  end
		  
		  if ~isreal(coh_data)
			error('Coherence matrix is not real')
		  end
					
		  coh_mat = coh_data;
		  
		  %
		  % Scale coherence matrix to have representation as int16
		  %
		  
		  coh_data = coh_data * 500;

		  iref_str = int2str(cf_list(iref));
		  if nwin_files == 1
			outfname1 = [new_root iref_str];
          		fid_mat = fopen(outfname1,'ab');
         		fseek(fid_mat, 0, 'eof');
		          num_written = fwrite(fid_mat,coh_data','int16');
			fclose(fid_mat);
		  else
	          for fref = 1:nwin_files
        		    text_ref = int2str(fref);
            			outfname1 = [new_root iref_str '_w' text_ref];
            			fid_mat = fopen(outfname1,'ab');
                		fseek(fid_mat, 0, 'eof');
				num_written = fwrite(fid_mat,coh_data','int16');
                  if num_written~= size(coh_data,1)*size(coh_data,2)
                    error('write failed')
                  end;

                		fclose(fid_mat);
          	  end
		  end
%		  if nwin_files > 1	
%		    num_written = fwrite(fid_mat(iref, iobs),coh_data','int16');
%		  else
%		    num_written = fwrite(fid_mat(iref, 1),coh_data','int16');
%		  end
		    
		  if num_written~= size(coh_data,1)*size(coh_data,2)
		    error('write failed')
		  end;
	    end;  %for iref =
	  end  % for iobs
	end; %for ifile =
end; % for icell
fclose('all');

format loose
status = 0;
return

