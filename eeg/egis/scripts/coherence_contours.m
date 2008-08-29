function status =coherence_contours(freqs,obs,array,coherencefname);
%status =coherence_contours(freqs,obs,array,coherencefname);
%
%generates contour map of coherence versus interelectrode distance and frequency
%freqs = frequencies to plot, e.g. [3:30] limits the plot to 3 to 30 Hz  
%obs = observations in the coherence file to plot (defaults to all)
%array = either an array name,e.g., 'all','hemisphere','quadrant','frontback',
%	or 'oned' or a list of channels
%coherencefname = coherencefilename (can be gui by skipping this param)
%

if nargin < 4
	[fid] = get_fid('rb');
else
	fid = fopen(coherencefname,'rb');
end;
version = fread(fid,1,'int16');
if version ~= -2
	error('this is not a coherence file');
end;
 [nfiles,obs_labels,Epoch,Window_Length,NEpoch,Nbad_chan,bad_chan,ref_flag,reference,NChan,NFreq] = rd_anal_hdr(fid);
frewind(fid);
if nargin == 0
	freqs = [0:40];
	obs = [1:nfiles];
	array = 'all';
	lines = 'black';
end;
if nargin == 1
	obs = [1:nfiles];
	array = 'all';
	lines = 'black';
end;		
if nargin == 2
	array = 'all';
	lines = 'black';
end;	

if isempty(array)
	array = 'all';
end
if isempty(obs)
	obs = [1:nfiles];
end;
if isempty(freqs)
	freqs = [0:40];
end;

ch_pair_indices;
coherencemax = 1;
iline = 1;
[xelec,yelec,zelec] = electrodes(129);
twod_dist = real(twod_pos(xelec,yelec,zelec,9.2));
if isstr(array)
	[arrayg,gnames] = arrays(array,[]);
else
	arrayg = array;
end
for g = 1:size(arrayg,1)
	for io = 1:size(obs,2)
		figure
		[npairs, pairs, pairnames] = look_up_pairs(array,g,bad_chan,obs(io));
		version = fread(fid,1,'int16');
		[nfiles,obs_labels,Epoch,Window_Length,NEpoch,Nbad_chan,bad_chan,ref_flag,reference,NChan,NFreq] = rd_anal_hdr(fid);
		coherence_array = real(get_eeg_coherence(fid,freqs,obs(io),nfiles,NChan,NFreq,Epoch));  		
		distance = [0:1:fix(max(twod_dist(pairs)))];                              
		coh_fit = zeros(size(freqs,2),size(distance,2));
		for j = 1:size(freqs,2)		

			p_coh = polyfit(twod_dist(pairs),coherence_array(j,pairs),4);
			coh_fit(j,:) = real(polyval(p_coh,distance));	
		end;		
		distance = [0:1:fix(max(twod_dist(pairs)))];
		hottie = hot;
		hotties = hottie(33:64,:);
		colormap(hotties);
		contourf(distance,freqs,coh_fit,[1:-0.1:0.1]),colorbar
		hold on, contour(distance,freqs,coh_fit,[0.1:0.1:0.2],'w--')
		hold on, contour(distance,freqs,coh_fit,[0.3:0.1:1],'k-')
		axis([0 max(twod_dist(pairs)) 0 max(freqs)]);
		ylabel('Frequency(Hz)')
		xlabel('Interelectrode Distance (cm)')
		title(['Coherence versus Distance and Frequency: '  deblank(setstr(obs_labels(obs(io),:))) ' ' pairnames]);
	end
end





























