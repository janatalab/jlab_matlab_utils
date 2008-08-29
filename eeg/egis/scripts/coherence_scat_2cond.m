function [chpair_coh] =coherence_scatter_2cond(freqs,obs,array,coh_crit,coherencefname);
% [chpair_coh] =coherence_scatter_2cond(freqs,obs,array,coh_crit,...
% coherencefname);
%freqs = frequencies to plot, e.g. [3 30] limits the plot at 3 and 30 Hz  
%obs = observations in the power file to plot (defaults to all)
%array = either an array name,e.g., 'all','hemisphere','quadrant','frontback',
%	or 'oned' or a list of channels
% Restrictions nobs <= 2. If nobs == 1 and multiple frequencies 
% are listed then freqs < = 2. If nobs > 1 obs go on same plot one for
% each freq.  if nobs == 1 all the freqs go on the same plot
% coh_crit = critical value for coherence change (** will default to 95%
%	confidence intervals but not implemented yet)
%coherencefname = coherencefilename (can be gui by skipping this param)
%
% chpair_coh = differnce data (see Manual for interpretation)
%
if nargin < 5
	[fid] = get_fid('rb');
else
	fid = fopen(coherencefname,'rb');
end;
version = fread(fid,1,'int16');
if version ~= -2
	error('this is not a coherence file');
end;
 [nfiles,obs_labels,Epoch,Window_Length,NEpoch,Nbad_chan,bad_chan,ref_flag,reference,NChan,NFreq] = rd_anal_hdr(fid);
if nargin == 0
	freqs = [0:40];
	obs = [1:nfiles];
	array = 'all';
	coh_crit = 0.1;
end;
if nargin == 1
	obs = [1:nfiles];
	array = 'all';
	coh_crit = 0.1;
end;		
if nargin == 2
	array = 'all';
	coh_crit = 0.1;
end;	
if nargin == 3
	coh_crit = 0.1;
end;
if isempty(coh_crit)
	coh_crit = 0.1;
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
if size(obs,2) > 2
	fclose('all')
	error('too many observations > 3 to be plotted in a single plot')
end;
if size(obs,2) == 1
	if size(freqs,2) > 2 
		error('too many frequencies for a single plot')
	end;
end;
lines = 'color';
ch_pair_indices;
if strcmp(lines,'black')
line_type = ['k- ';'k--';'k-.'];
scat_type = ['ko';'k.';'ko'];
else	
line_type = ['k-';'r-';'b-'];
scat_type = ['ko';'r.';'bo'];
end;
coherencemax = 1;
iline = 1;
[xelec,yelec,zelec] = electrodes(129);
twod_dist = twod_pos(xelec,yelec,zelec,9.2);
coherence_array = real(get_eeg_coherence(fid,freqs,obs,nfiles,NChan,NFreq,Epoch));
[dist,bcohpotential,bcohlaplacian] = baseline_coherence(0.2,80,8,8.2,8.7,9.2,7.8);
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
				for io = 1:size(obs,2)
					if ~(ref_flag(obs(io)) == 5|ref_flag(obs(io)) == 7)
						hold on,bar(dist,bcohpotential,'w-');
					elseif ref_flag(obs(io)) == 5 
						hold on,bar(dist,bcohlaplacian,'w-');
					end;
				end
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
%				if j == 1
%				hold on, plot(real(twod_dist(pairs)),coherence_array((j-1)*size(freqs,2)+i,pairs),scat_type(2,:))
%				end;
				if j == 2
				hold on, plot(real(twod_dist(pairs)),coherence_array((j-1)*size(freqs,2)+i,pairs),scat_type(iline,:))
				diff_pairs = find(abs(coherence_array(size(freqs,2)+i,:) - coherence_array(i,:)) > coh_crit);
				mask = zeros(1,size(coherence_array,2));
				mask(pairs) = ones(1,size(pairs,2));
				mask_diff = zeros(1,size(coherence_array,2));
				mask_diff(diff_pairs) = ones(1,size(diff_pairs,2));
				new_mask = mask_diff.*mask;
				pairs_diff = find(new_mask);
				hold on, plot(real(twod_dist(pairs_diff)),coherence_array((j-1)*size(freqs,2)+i,pairs_diff),scat_type(3,:))
				end
%				p_coh = polyfit(twod_dist(pairs),coherence_array((j-1)*size(freqs,2)+i,pairs),4);
%				coh_fit = real(polyval(p_coh,[0:1:fix(max(twod_dist(pairs)))]));
%				hold on, plot([0:1:max(twod_dist(pairs))],coh_fit,line_type(iline,:)), 
				hold on, axis([0 max(twod_dist(pairs))+1 0 1])
				title(['Coherence versus Distance:  @' int2str(freqs(i)) ' Hz ' hnames(g,:)]);
				text([max(twod_dist(pairs))-2.5],[1-0.1*j],[scat_type(iline+1,2) ' ' deblank(setstr(obs_labels(obs(j),:)))])
				xlabel('Interelectrode Distance (cm)')
				ylabel('Coherence')
 				iline = iline+ 1;
			end
		chpair_list = chpairindex_2chpair(pairs_diff);
		chpair_coh((g-1)*size(chpair_list,1)+1:g*size(chpair_list,1),:) = [g*ones(size(chpair_list,1),1) freqs(i)*ones(size(chpair_list,1),1) obs(1)*ones(size(chpair_list,1),1) obs(2)*ones(size(chpair_list,1),1) chpair_list twod_dist(pairs_diff)' coherence_array(i,pairs_diff)' coherence_array(i+size(freqs,2),pairs_diff)'];
		
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
		for io = 1:size(obs,2)
			if ~(ref_flag(obs(io)) == 5|ref_flag(obs(io)) == 7)
				hold on,bar(dist,bcohpotential,'w-');
			elseif ref_flag(obs(io)) == 5
				hold on,bar(dist,bcohlaplacian,'w-');
			end;
		end
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
%			chpair(pairels(1),pairels(2))
			pairnames = hnames(g,:);
%			if j == 1
%			hold on, plot(real(twod_dist(pairs)),coherence_array(j,pairs),scat_type(iline,:))
%			end;
			if j == 2
			hold on, plot(real(twod_dist(pairs)),coherence_array(j,pairs),scat_type(2,:))
			diff_pairs = find(abs(coherence_array(2,:) - coherence_array(1,:)) > coh_crit);
			mask = zeros(1,size(coherence_array,2));
			mask(pairs) = ones(1,size(pairs,2));
			mask_diff = zeros(1,size(coherence_array,2));
			mask_diff(diff_pairs) = ones(1,size(diff_pairs,2));
			new_mask = mask_diff.*mask;
			pairs_diff = find(new_mask);
			hold on, plot(real(twod_dist(pairs_diff)),coherence_array(j,pairs_diff),scat_type(3,:))
			end
%			p_coh = polyfit(real(twod_dist(pairs)),coherence_array(j,pairs),4);
%			coh_fit = real(polyval(p_coh,[0:1:fix(max(twod_dist(pairs)))]));
%			hold on, plot([0:1:max(twod_dist(pairs))],coh_fit,line_type(iline,:)), 
			hold on, axis([0 max(twod_dist(pairs))+1 0 1])
			title(['Coherence versus Distance: ' deblank(setstr(obs_labels(obs(1),:))) ' ' hnames(g,:)]);
			text([max(twod_dist(pairs))-2.5],[1-0.1*j],[scat_type(iline+1,2) ' ' int2str(freqs(j)) ' Hz'])
			xlabel('Interelectrode Distance (cm)')
			ylabel('Coherence')
			iline = iline+ 1;
		end
		chpair_list = chpairindex_2chpair(pairs_diff);
		chpair_coh((g-1)*size(chpair_list,1)+1:g*size(chpair_list,1),:) = [g*ones(size(chpair_list,1),1) obs(1)*ones(size(chpair_list,1),1) freqs(1)*ones(size(chpair_list,1),1) freqs(2)*ones(size(chpair_list,1),1) chpair_list twod_dist(pairs_diff)' coherence_array(1,pairs_diff)' coherence_array(2,pairs_diff)'];
		
	end
end
status = 1;





























