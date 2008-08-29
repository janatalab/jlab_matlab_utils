function [npairs, pairs, pairnames] = look_up_pairs(array,iarray,bad_chan,obs);
ch_pair_indices;
if isstr(array)
	[arrayh,hnames] = arrays(array,bad_chan(obs,find(bad_chan(obs,:))));
else
	arrayh = array;
end
pairels = arrayh(iarray,find(arrayh(iarray,:)));
npairs = size(pairels,2)*(size(pairels,2)+1)/2;
pairs = zeros(1,npairs);
ipairs = 1;
for ip = 1:size(pairels,2)
	for jp=ip:size(pairels,2)
		pairs(ipairs) = chpair(pairels(ip),pairels(jp));
		ipairs = ipairs + 1;
	end;
end;
pairnames = hnames(iarray,:);
