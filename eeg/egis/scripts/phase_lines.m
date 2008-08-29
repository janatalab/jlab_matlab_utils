function status =phase_lines(freqs,obs,array,lines,phase_type,phasefname);
%status =phase_lines(freqs,obs,array,lines,phasefname);
%
%freqs = frequencies to plot, e.g. [3 30] limits the plot at 3 and 30 Hz  
%obs = observations in the power file to plot (defaults to all)
%array = either an array name,e.g., 'all','hemisphere','quadrant','frontback',
%	or 'oned' or a list of channels
% Restrictions nobs <= 3. If nobs == 1 and multiple frequencies 
% are listed then freqs < = 3. If nobs > 1 obs go on same plot one for
% each freq.  if nobs == 1 all the freqs go on the same plot
%lines = 'black' or 'color'
%phasefname = phasefilename (can be gui by skipping this param)
%
%note also draws in solid black phase due to uncorrelated sources
%in the four spheres model.  If 'laplacian' is a reference the curve
%for laplacians is drawn as well. 
if nargin < 6
	[fid] = get_fid('rb');
else
	fid = fopen(phasefname,'rb');
end;
version = fread(fid,1,'int16');
if version ~= -2
	error('this is not a phase file');
end;
 [nfiles,obs_labels,Epoch,Window_Length,NEpoch,Nbad_chan,bad_chan,ref_flag,reference,NChan,NFreq] = rd_anal_hdr(fid);
if nargin == 0
	freqs = [0:40];
	obs = [1:nfiles];
	array = 'all';
	lines = 'black';
	phase_type = 'angles';
end;
if nargin == 1
	obs = [1:nfiles];
	array = 'all';
	lines = 'black';
	phase_type = 'angles';
end;		
if nargin == 2
	array = 'all';
	lines = 'black';
	phase_type = 'angles';
end;	
if nargin == 3
	lines = 'black';
	phase_type = 'angles';
end;
if nargin == 4
	phase_type = 'angles';
end;
if isempty(lines)
	lines = 'black';
end
if isempty(array)
	array = 'all';
end
if isempty(obs)
	obs = [1:nfiles];
end;
if isempty(freqs)
	freqs = [0:40];
end;
if isempty(phase_type)
	phase_type = 'angles';
end;
if size(obs,2) > 3
	fclose('all')
	error('too many observations > 3 to be plotted in a single plot')
end;
if size(obs,2) == 1
	if size(freqs,2) > 3 
		error('too many frequencies for a single plot')
	end;
end;
ch_pair_indices;
if strcmp(lines,'black')
line_type = ['k- ';'k--';'k-.'];
scat_type = ['kx';'ko';'k.'];
else	
line_type = ['b-';'r-';'k-'];
scat_type = ['kx';'ko';'k.'];
end;

iline = 1;
[xelec,yelec,zelec] = electrodes(129);
twod_dist = twod_pos(xelec,yelec,zelec,9.2);
phase_array = real(get_eeg_coherence(fid,freqs,obs,nfiles,NChan,NFreq,Epoch));
if strcmp(phase_type,'cos')
	phase_array = cos(phase_array);
	phasemax = 1;
else
	phasemax = pi;
end;
if isstr(array)
	[arrayg,gnames] = arrays(array,[]);
else
	arrayg = array;
end
if size(obs,2) > 1
	for g = 1:size(arrayg,1)
			for i= 1:size(freqs,2)
				figure
				hold on
				
				iline = 1;
				
				for j = 1:size(obs,2)
				
					if isstr(array)
						[arrayh,hnames] = arrays(array,bad_chan(obs(j),find(bad_chan(obs(j),:))));
					else
						arrayh = array;
					end
				pairels = arrayh(g,find(arrayh(g,:)));
				npairs = size(pairels,2)*(size(pairels,2)+1)/2;
				pairs = zeros(1,npairs);
				ipairs = 1;
				for ip = 1:size(pairels,2)
					for jp=ip:size(pairels,2)
						pairs(ipairs) = chpair(pairels(ip),pairels(jp));
						ipairs = ipairs + 1;
					end;
				end;
				pairnames = hnames(g,:);
%				hold on, plot(real(twod_dist(pairs)),phase_array((j-1)*size(freqs,2)+i,pairs),scat_type(iline,:))
				p_coh = polyfit(twod_dist(pairs),phase_array((j-1)*size(freqs,2)+i,pairs),4);
				coh_fit = real(polyval(p_coh,[0:1:fix(max(twod_dist(pairs)))]));
				hold on, plot([0:1:max(twod_dist(pairs))],coh_fit,line_type(iline,:)), axis([0 max(twod_dist(pairs))+1 -phasemax*1.1 phasemax*1.1])
				title(['Phase versus Distance:  @' int2str(freqs(i)) ' Hz ' hnames(g,:)]);
				hold on, plot([max(twod_dist(pairs))-3  max(twod_dist(pairs))-2 ],[1-0.1*j 1-0.1*j],line_type(iline,:))
				text([max(twod_dist(pairs))-2.5],[1-0.1*j],[deblank(setstr(obs_labels(obs(j),:)))])
				xlabel('Interelectrode Distance (cm)')
				ylabel('Phase')
 				iline = iline+ 1;
			end
		end
	end
else
	if isstr(array)
		[arrayh,hnames] = arrays(array,bad_chan(obs,find(bad_chan(obs,:))));
	else
		arrayh = array;
	end
	for g = 1:size(arrayh,1)
		iline = 1;
		figure
		hold on
		
		for j = 1:size(freqs,2)
			pairels = arrayh(g,find(arrayh(g,:)));
			
			npairs = size(pairels,2)*(size(pairels,2)+1)/2;
			pairs = zeros(1,npairs);
			ipairs = 1;
			for ip = 1:size(pairels,2)
				for jp=ip:size(pairels,2)
					pairs(ipairs) = chpair(pairels(ip),pairels(jp));
					ipairs = ipairs + 1;
				end;
			end;
			chpair(pairels(1),pairels(2))
			pairnames = hnames(g,:);
			
%			hold on, plot(real(twod_dist(pairs)),phase_array(j,pairs),scat_type(iline,:))
			p_coh = polyfit(real(twod_dist(pairs)),phase_array(j,pairs),4);
			coh_fit = real(polyval(p_coh,[0:1:fix(max(twod_dist(pairs)))]));
			hold on, plot([0:1:max(twod_dist(pairs))],coh_fit,line_type(iline,:)), axis([0 max(twod_dist(pairs))+1 -phasemax*1.1 phasemax*1.1])
			title(['Phase versus Distance: ' deblank(setstr(obs_labels(obs(1),:))) ' ' hnames(g,:)]);
			hold on, plot([max(twod_dist(pairs))-3 max(twod_dist(pairs))-2],[1-0.1*j 1-0.1*j],linetype(iline,:))
			text([max(twod_dist(pairs))-2.5],[1-0.1*j],[int2str(freqs(j)) ' Hz'])
			xlabel('Interelectrode Distance (cm)')
			ylabel('Phase')
			iline = iline+ 1;
		end
	end
end
status = 1;





























