function status = make_factor(csdmfiles,outfname,win,freqlist,nfactor,bad_chan,output);
%status = make_factor(csdmfiles,outfname,win,binmin,binmax,nfactor,bad_chan,output);
%
%	this function makes factor files derived files from csdm files
%	the names of the csdm files must be stored in an array where each row 
%		is a filename - csdmfiles.   
%	outfname = root name of output file factor # is appended
%	win = which windows (observation number) to use in each file 
%	COMING SOON	ifwin = -1 generate one file for each observation
%			if win = 0 generate one file fo reach subject
%	freqlist = frequency list to analyze
%	nfactor = #of factor files to create
%	bad_chan = badchannel list 
%	output = 'cos' - cosine of phase 'ang' phase angle defaults to 'cos'
%

if ~(nargin==6|nargin==7|nargin == 5)
	error('improper arguments')
end;

if nargin == 5
	bad_chan = [];
	output = 'cos';
elseif nargin == 6
	output = 'cos';
end;


%
% read in a csdm header in order to make an output file header
%
ave_hdr_offsets_v;
%csdm_fid = fopen(csdmfiles(1,:),'r');
%[csdm_fhdr,csdm_chdr,ename,czeros,cgains,cnames,fcom,ftext, coff]=rd_csdm_hdr_v(csdm_fid);
%fclose(csdm_fid);
%[cell_data_offsets]=get_csdm_cell_off(csdm_fhdr, csdm_chdr);%
%
%if win > csdm_chdr(1,NObs)
%	error('window out of range')
%end;
%
%ave_fhdr = csdm_fhdr;
% Copy average file subject specifics array
%ave_chdr = concat_csdm_chdr(csdmfiles, win);

%if csdm_fhdr(NChan) == 8385
%	ave_fhdr(NChan) = 129;
%elseif csdm_fhdr(NChan) == 2145
%	ave_fhdr(NChan) = 65;
%elseif csdm_fhdr(NChan) == 33153;
%	ave_fhdr(NChan) = 257;
%else
%	error('Nchan is unacceptable in this csdm file')
%end;
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
        ave_chdr(c,NPoints) = size(freqlist,2);
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
for i = 1:nfactor
	for f = 1:nwin_files
	  if nwin_files > 1 & win ~= []
    		outfname1 = [out_root '_mag_' int2str(win(f)) '_ft_' int2str(i)];
    		outfname2 = [outfname '_' output '_' int2str(win(f)) '_ft_' int2str(i)];
  	  else
	    outfname1 = [out_root '_mag_ft_' int2str(i)];
    	    outfname2 = [out_root '_' output '_ft_' int2str(i)];
       	  end
	
  	  ave_fid = fopen(outfname1,'wb');
  	  ave_fid2 = fopen(outfname2,'wb');
  	  wt_ave_hdr_v(ave_fid,ave_fhdr,ave_chdr,ave_ename,ave_cnames,ave_fcom,ave_ftext);
  	  wt_ave_hdr_v(ave_fid2,ave_fhdr,ave_chdr,ave_ename,ave_cnames,ave_fcom,ave_ftext);

  	fclose(ave_fid);
  	fclose(ave_fid2);
	end
end;

ch_pair_indices;

[cell_data_offsets]=get_csdm_cell_off(csdm_fhdr, csdm_chdr);

file_mask = ones(size(csdmfiles,1),ave_fhdr(NChan));
if (nargin == 7|nargin == 8)
	for i=1:size(bad_chan,1)
		bchan = find(bad_chan(i,:));
		if bchan ~= []
			file_mask(i,bad_chan(bchan)) = zeros(1,size(bchan,2));
		end;
	end;
end;
%for i=1:nfactor
%	outfname1 = [outfname '_mag_' int2str(i)]
%	outfname2 = [outfname '_' output '_' int2str(i)]
%	avefid(i) = fopen(outfname1,'wb');
%	avefid(i+nfactor) = fopen(outfname2, 'wb');
%	% Write out the average file header
%	wt_ave_hdr_v(avefid(i),ave_fhdr,ave_chdr,ave_ename,ave_cnames,ave_fcom,ave_ftext);
%	fseek(avefid(i), 0, 'eof');
%        wt_ave_hdr_v(avefid(i+nfactor),ave_fhdr,ave_chdr,ave_ename,ave_cnames,ave_fcom,ave_ftext);
%        fseek(avefid(i+nfactor), 0, 'eof');
%	fclose(avefid(i));
%	fclose(avefid(i+nfactor));
%end;

ch_pair_indices;
outfname3 = [outfname '_eig'];
eig_fid = fopen(outfname3,'wb');
fwrite(eig_fid,[csdm_fhdr(NCells) size(csdmfiles,1) ave_chdr(1,NPoints) nwin],'int16');
for icell = 1:csdm_fhdr(NCells)
	for ifile = 1:size(csdmfiles,1)
	   for f = 1:nwin
		[csdm_data, csdm_fhdr, csdm_chdr] = read_csdm(csdmfiles(ifile,:), icell, win(f));
			
		csdm_mat = zeros(ave_fhdr(NChan));
		for jfreq = 1:size(freqlist,2)
			ifreq= freqlist(jfreq);
			for ichan =1:ave_fhdr(NChan)
				for jchan = ichan:ave_fhdr(NChan);
					csdm_mat(ichan,jchan) = csdm_data(ifreq,chpair(ichan,jchan));
					csdm_mat(jchan,ichan) = conj(csdm_mat(ichan,jchan));
				end;
			end;
			good_chan = find(file_mask(ifile,:));
			csd = zeros(size(good_chan,2));
			csd = csdm_mat(good_chan,good_chan);
		
			[Vec_csd, val_csd] = eig(csd,'nobalance');
%			val_csd([1:3 127:129],[1:3 127:129])
			if ~isreal(val_csd)
%				disp('eigenvalues are not real, fixing')
%				save val_csdm
%				val_csd(size(good_chan,2),size(good_chan,2))
				val_csd = real(val_csd);
			end;
			vec_csd = vec_csd';
			abs_vec_csd = 500*abs(vec_csd);
			fwrite(eig_fid,[icell ifreq size(val_csd,1) f],'int16');
%			fwrite(eig_fid,ifile,'int16');
%			fwrite(eig_fid,ifreq,'int16');
%			fwrite(eig_fid,size(val_csd,1),'i
				fwrite(eig_fid,diag(val_csd),'float');
		
			phase_data = vec_csd;
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

			if ~isreal(phase_data)
                  		error('Phase matrix is not real')
                	end
%			if nwin_files > 1 & win ~= []
%			    outfname1 = [out_root '_mag_' int2str(win(f))];
%    			    outfname2 = [outfname '_' output '_' int2str(win(f))];
%  			else
%    			   outfname1 = [out_root '_mag'];
% 			   outfname2 = [out_root '_' output];
%  			end


			for ifac = 1:nfactor
				mag = zeros(1, ave_fhdr(NChan));
				ph = zeros(1,ave_fhdr(NChan));
				mag(1,good_chan) = abs_vec_csd(ifac,:);
				ph(1,good_chan) = phase_data(ifac,:);
%			        outfname1 = [outfname '_mag_' int2str(ifac)];
%        			outfname2 = [outfname '_' output '_' int2str(ifac)];
	                        if nwin_files > 1 & win ~= []
        	                    outfname1 = [out_root '_mag_' int2str(win(f)) '_ft_' int2str(ifac)];
                	            outfname2 = [outfname '_' output '_' int2str(win(f)) '_ft_' int2str(ifac)];
                        	else
                           		outfname1 = [out_root '_mag_ft_' int2str(ifac)];
                           		outfname2 = [out_root '_' output '_ft_' int2str(ifac)];
                        	end

        			avefid1 = fopen(outfname1,'ab');
        			avefid2= fopen(outfname2, 'ab');

				num_written = fwrite(avefid1, mag','int16');
				if num_written~= size(mag,2)
					error('write out croaked')
				end;
				num_written = fwrite(avefid2, ph','int16');
                        	if num_written~= size(ph,2)
                                	error('write out croaked')
                        	end;
				fclose(avefid1);
				fclose(avefid2);
			end	
		end;
	  end;
	end; %for ifile =
end; % for icell
fclose('all');



