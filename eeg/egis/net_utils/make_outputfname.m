function powerfname = make_outputfname(csdmfname,output_type)

for i = 2:size(csdmfname,1)
	subj_check(i) = strcmp(csdmfname(1,1:3),csdmfname(i,1:3));
	group_check(i) =  strcmp(csdmfname(1,5:5),csdmfname(i,5:5));
	cnd1(i) = strcmp(csdmfname(1,7:7),csdmfname(i,7:7));
	cnd2(i) = strcmp(csdmfname(1,8:8),csdmfname(i,8:8));
	cnd3(i) = strcmp(csdmfname(1,9:9),csdmfname(i,9:9));
	cnd4(i) = strcmp(csdmfname(1,10:10),csdmfname(i,10:10));
	cnd5(i) = strcmp(csdmfname(1,11:11),csdmfname(i,11:11));
	ref_check(i) = strcmp(csdmfname(1,13:16),csdmfname(i,13:16));
end;
iff = 1;
if sum(subj_check) == size(csdmfname,1) - 1
	 powerfname(1,iff:iff+3) = [csdmfname(1,1:3) '_'];
	iff = iff +4;
end
if sum(group_check) == size(csdmfname,1) - 1
	powerfname(1,iff:iff+1) = [csdmfname(1,5) '_'];
	iff = iff+2;
end;
if (sum(cnd1) == size(csdmfname,1) - 1) & (csdmfname(1,7) ~= 0)
	powerfname(1,iff) = [csdmfname(1,7)];	
	iff = iff+1;
end;
if (sum(cnd2) == size(csdmfname,1) - 1) & (csdmfname(1,8) ~= 0)
	powerfname(1,iff) = [csdmfname(1,8)];
	iff = iff + 1;
end
if (sum(cnd3) == size(csdmfname,1) - 1) & (csdmfname(1,9) ~= 0)
	powerfname(1,iff) = [csdmfname(1,9)];
	iff = iff + 1;
end
if (sum(cnd4) == size(csdmfname,1) - 1) & (csdmfname(1,10) ~= 0)
	powerfname(1,iff) = [csdmfname(1,10)];
	iff = iff + 1;
end
if (sum(cnd5) == size(csdmfname,1) - 1) & (csdmfname(1,11) ~= 0)
	powerfname(1,iff) = [csdmfname(1,11)];
	iff = iff + 1;
end
if sum(ref_check) == size(csdmfname,1) - 1
	powerfname(1,iff:iff+4) = ['.' csdmfname(1,13:16)];
end;
if strcmp(output_type,'power')
	powerfname = [powerfname '.pow'];
elseif strcmp(output_type,'phase')
	powerfname = [powerfname '.pha'];
elseif	strcmp(output_type,'coherence')
	powerfname = [powerfname '.coh'];
else
	disp('unknown type no type flag appended')
end;


