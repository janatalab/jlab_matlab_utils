function chpair_list = chpairindex_2chpair(pairindex);
ch_pair_indices;
for i =1:size(pairindex,2)
	[chpair_list1, chpair_list2] = find(chpair == pairindex(i));
	chpair_list(i,1:2) = [chpair_list1(1) chpair_list2(1)]; 
end;
%icount = 1;
%for i = 1:size(chpair_list1,1)
%	if chpair_list1(i) < chpair_list2(i)
%		chpair_list(icount,1:2) = [chpair_list1(icount) chpair_list2(icount)];
%		icount = icount + 1;
%	end;
%end;
%size(chpair_list,1)
%size(pairindex,2)
if size(chpair_list,1) ~= size(pairindex,2)
	error('houston we have a problem')
end;

