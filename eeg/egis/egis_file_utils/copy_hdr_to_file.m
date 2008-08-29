function status=copy_hdr_to_file(infname,outfname)
%status=copy_hdr_to_file(infname,outfname)
%
%Writes and EGIS file header to a new file.

%Use this code for batch processing
infid=fopen(infname);
if infid == -1, error('Fuck up opening input file.');, end

[fhdr,chdr,ename,czeros,cgains,cnames,fcom,ftext]=rd_ses_hdr_v(infid);

%Use this code for batch processing
outfid=fopen(outfname,'w');
if outfid == -1
	if infid ~= -1 
		status=fclose(infid)
	end
	error('Fuck up opening output file.');
end

status=wt_ses_hdr_v(outfid,fhdr,chdr,ename,czeros,cgains,cnames,fcom,ftext);

status=fclose('all');