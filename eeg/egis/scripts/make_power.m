function power_data = make_power(csdmfiles,outfname,win,output);
%power_data = make_power(csdmfiles,outfname,win,output);
%
%	this function makes power files derived files from csdm files
%	the names of the csdm files must be stored in an array where each row 
%		is a filename - csdmfiles.   
%	outfname = name of output file
%	win = 	which window (observation number) to use in each file 
%			if win == [], all observations are used
%			if length(win) > 0 and there is more than one csdmfile specified, a
%           separate file is created for each window, with NObs in each file
%			corresponding to the number of input csdmfiles
%	output = 'amp' for amplitude, 'pow' for power, and 'log' for log
%		defaults to 'amp'
% 

% Original version by RS
% Modification history:
%  Introduced get_csdm_hdr and csdmhdr2avehdr -- PJ
%
%  02/05/96 PJ added 'norm' option which returns proportion variance
%              explained scaled by 500.  This eliminates the risk of int16 overflow
%              as when using 'pow'.

if ~(nargin==3|nargin==4)
	error('improper arguments')
end;

if nargin == 3
	output ='amp';
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
	ave_chdr(c,NPoints) = csdm_chdr(c,NPoints)/2;
end

%
% Determine the number of output files
%

if size(csdmfiles,1) > 1 & length(win) > 1
  nwin_files = length(win);
else
  nwin_files = 1;
end

nwin = length(win);

%
% Write out the average file header -- same header if multiple files
%

out_root = outfname;
for f = 1:nwin_files
  if nwin_files > 1 & win ~= []
    outfname1 = [out_root '_' int2str(win(f))];
  else
    outfname1 = out_root;
  end
end
  
  ave_fid = fopen(outfname1,'wb');
  wt_ave_hdr_v(ave_fid,ave_fhdr,ave_chdr,ave_ename,ave_cnames,ave_fcom,ave_ftext);
  fclose(ave_fid);
end
if ave_fhdr(NChan) == 129
ch_pair_indices;
else
ch_pair_indices_65;
end;

[cell_data_offsets]=get_csdm_cell_off(csdm_fhdr, csdm_chdr);

for icell = 1:csdm_fhdr(NCells)
  disp(['Working on cell ' int2str(icell)]);
	for ifile = 1:size(csdmfiles,1)
	  disp(['  Reading from file: ' csdmfiles(ifile,:)]);
	  for iwin = 1:nwin
	    disp(['    Window ' int2str(iwin) '  (' int2str(win(iwin)) ')']);
		[csdm_data, csdm_fhdr, csdm_chdr] = read_csdm(csdmfiles(ifile,:), icell, win(iwin));
		npts = ave_chdr(icell,NPoints);
		power_data = zeros(ave_chdr(icell,NPoints),ave_fhdr(NChan));
		for ichan = 1:ave_fhdr(NChan)
			power_data(:,ichan) = csdm_data(:,chpair(ichan,ichan));
		end;
		
		% Reality check
		testreal = isreal(sum(power_data));
		if testreal ~= 1
			power_data(11,73)
			error('complex numbers')
		end;
		
		if strcmp(output,'amp')
		 	power_data = 500*sqrt(power_data);
		elseif strcmp(output,'pow')	
			power_data = 500*power_data;
		  while any(max(max(power_data)) > 2^15)
		    disp('Warning: Int16 overflow -- scaling data down by factor of 500');
		    power_data = power_data / 500;
          end
		elseif strcmp(output,'norm')
%		  power_data(2:npts,:) = power_data(2:npts,:) ./ (ones(npts-1,1) * sum(power_data(2:npts,:))) * 500;
		  power_data = power_data ./ (ones(npts,1) * sum(power_data)) * 500;
		else  % default to log10
			power_data = 500*log10(power_data);
		end

		if nwin_files > 1
		  outfname1 = [out_root '_' int2str(win(iwin))];
                else
		  outfname1 = out_root;
		end

		ave_fid = fopen(outfname1, 'ab');
		fseek(ave_fid, 0, 'eof');

		num_written = fwrite(ave_fid,power_data','int16');
		if num_written~= size(power_data,1)*size(power_data,2)
			error('write failed')
		end;

		fclose(ave_fid);
		
	  end; % for iwin = 
	end; %for ifile =
end; % for icell
fclose('all');


