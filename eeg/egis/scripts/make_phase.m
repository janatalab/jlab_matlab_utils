function status = make_phase(csdmfiles,outfname,win,rel_chan,output);
%status = make_phase(csdmfiles,outfname,win,rel_chan,output);
%
%	this function makes phase files derived files from csdm files
%	the names of the csdm files must be stored in an array where each row 
%		is a filename - csdmfiles.   
%	outfname = name of output file with rel_chan appended
%	win = which window (observation number) to use in each file 

%	COMING SOON	if win = -1 generate one file for each observation
%			if win = 0 generate one file fo reach subject
%	rel_chan = index channels for phase
%	output = 'ang' for angle, 'cos' for cosine phase
%		defaults to 'cos'
%

if ~(nargin==5|nargin==4)
	error('improper arguments')
end;

if nargin == 4
	output ='cos';
end
%
% read in a csdm header in order to make an output file header
%
ave_hdr_offsets_v;
%csdm_fid = fopen(csdmfiles(1,:),'r');
%[csdm_fhdr,csdm_chdr,ename,czeros,cgains,cnames,fcom,ftext, coff]=rd_csdm_hdr_v(csdm_fid);
%
%[cell_data_offsets]=get_csdm_cell_off(csdm_fhdr, csdm_chdr);
%
%if win > csdm_chdr(1,NObs)
%	error('window out of range')
%end;%
%
%ave_fhdr = csdm_fhdr;
%if csdm_fhdr(NChan) == 8385
%	ave_fhdr(NChan) = 129;
%elseif csdm_fhdr(NChan) == 2145
%	ave_fhdr(NChan) = 65;
%elseif csdm_fhdr(NChan) == 33153;
%	ave_fhdr(NChan) = 257;
%else
%	error('Nchan is unacceptable in this csdm file')
%end;
 

% Copy average file subject specifics array
%ave_chdr = concat_csdm_chdr(csdmfiles, win);
%
%ave_ename = ename; ave_fcom = fcom; ave_ftext = ftext;
%ave_fhdr(LastDone) = size(csdmfiles,1);
%ave_cnames = cnames;

% Copy the experiment info from the session file to the average file tspecs
%[ave_chdr] = sesfhdr2avetspecs(ses_fhdr, ave_chdr);

%Fix header 
%for c = 1:csdm_fhdr(NCells)
%        ave_chdr(c,NObs) = size(csdmfiles,1);
%	ave_chdr(c,NPoints) = csdm_chdr(c,NPoints)/2;
%end
%for i = 1:size(rel_chan,2)
%	outfname1 = [outfname '_' int2str(rel_chan(i))];
%	ave_fid(i) = fopen(outfname1,'wb');
%	% Write out the average file header
%	wt_ave_hdr_v(ave_fid(i),ave_fhdr,ave_chdr,ave_ename,ave_cnames,ave_fcom,ave_ftext);
%%	fseek(ave_fid(i), 0, 'eof');
%end;
%fclose('all');

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
  nwin_files = length(win)*size(rel_chan,2);
else
  nwin_files = size(rel_chan,2);
end;

nwin = length(win);

out_root = [outfname '_' output];
if size(csdmfiles,1) > 1 & length(win) > 1
	for f = 1:size(win,2)
		for r=1:size(rel_chan,2)
		    outfname1 = [out_root '_w_' int2str(win(f)) '_r_' int2str(rel_chan(r))];
		      ave_fid = fopen(outfname1,'wb');
		     wt_ave_hdr_v(ave_fid,ave_fhdr,ave_chdr,ave_ename,ave_cnames,ave_fcom,ave_ftext);
  		     fclose(ave_fid);

		end;
	end;

  else
	for f = 1:size(rel_chan,2)
  		outfname1 = [out_root '_r_' int2str(rel_chan(f))];
		 ave_fid = fopen(outfname1,'wb');
                     wt_ave_hdr_v(ave_fid,ave_fhdr,ave_chdr,ave_ename,ave_cnames,ave_fcom,ave_ftext);
                     fclose(ave_fid);
	end;

    
end






ch_pair_indices;
for icell = 1:csdm_fhdr(NCells)
	for ifile = 1:size(csdmfiles,1)
	   for iwin = 1:nwin
		[csdm_data, csdm_fhdr, csdm_chdr] = read_csdm(csdmfiles(ifile,:), icell, win(iwin));
		for irel = 1:size(rel_chan,2)
			phase_data = zeros(ave_chdr(icell,NPoints),ave_fhdr(NChan));
			for ichan = 1:ave_fhdr(NChan)
				phase_data(:,ichan) = csdm_data(1:ave_chdr(icell,NPoints),chpair(ichan,rel_chan(irel)));
				if ichan > rel_chan(irel)
					phase_data(:,ichan) = conj(phase_data(:,ichan));
				end
			end

			if output == 'cos'
		 		phase_data = 500*real(phase_data)./abs(phase_data);
			else
				real_sign = sign(real(phase_data));
				phase_data = angle(phase_data);
				ang_sign = sign(phase_data);
				for iph = 1:size(phase_data,1)
					for jph = size(phase_data,2)
						if real_sign(iph,jph) < 0
							phase_data(iph,jph) = phase_data(iph,jph) + pi*real_sign(iph,jph)*ang_sign(iph,jph);
						end;
					end;
				end;
				
		
				phase_data = 500*phase_data;
			end
 			phase_data(find(isnan(phase_data))) = zeros(size(find(isnan(phase_data))));
			if size(csdmfiles,1) > 1 & length(win) > 1
        			f = iwin;
                		r = irel;
                    		outfname1 = [out_root '_w_' int2str(win(f)) '_r_' int2str(rel_chan(r))];
                      		fid1 = fopen(outfname1,'ab');
                     		num_written = fwrite(fid1,phase_data','int16');
                        	fclose(fid1);
                        	if num_written~= size(phase_data,1)*size(phase_data,2)
                                	error('write failed')
                        	end;

                     	else

		        	f = irel;
                		outfname1 = [out_root '_r_' int2str(rel_chan(f))];
                 		fid1 = fopen(outfname1,'ab');
                     		num_written = fwrite(fid1,phase_data','int16');
                        	fclose(fid1);
                        	if num_written~= size(phase_data,1)*size(phase_data,2)
                                	error('write failed')
                        	end;

        		end;




%			outfname1 = [outfname '_' int2str(rel_chan(irel))];
%			fid1 = fopen(outfname1,'ab');
%			num_written = fwrite(fid1,phase_data','int16');
%			fclose(fid1);
%			if num_written~= size(phase_data,1)*size(phase_data,2)
%				error('write failed')			end;
		end
           end
	end; %for ifile =
end; % for icell
fclose('all');



