function status=make_gave_diff(cell_map, infname, outfname);
%function status=make_gave_diff(cell_map, infname, outfname);
%
%This function creates a difference file with name outfname
%from file infname using cell_map as a guide. If infname &
%outfname are omitted, this function interactively prompts
%for input and output files. The argument cell_map should
%contain one row for each difference you wish to calculate
%and put in the output file. The first column is the source
%cell from the input file, and the second column is the
%cell to be subtracted.
%

%Modification History
%
%Created Dec. 1, 1995 by Brian C. Rakitin
%

status = -1;

%Check number of input arguments
if ~((nargin == 1)|(nargin == 3))
	error('Number of input arguments must be either 1, or 3.');
end

%Initialize fids
destfid = -1;
srcfid = -1;
%First try batch mode
if nargin == 3
	srcfid= fopen(infname(i,:), 'r');
	if srcfid== -1
		error(['Could not open input file' infname '.']);
	end
	destfid = fopen(outfname, 'wb');
	if destfid == -1
		temp =fclose(srcfid);
	error(['Could not open output file ' outfname '.']);
	end
%Otherwise run interactive mode
elseif nargin <= 3
	while srcfid == -1
		[srcfid, infname]=get_fid('r','*.gave*', 'Open Average File:');
	end
	while destfid == -1
		[destfid, outfname]=put_fid('wb','new.gave_diff','Save GAVE Diff. File As:');
	end
end

ave_hdr_offsets_v;
[fhdr,chdr,ename,czeros,cgains,cnames,fcom,ftext, coff]=rd_egis_hdr_v(srcfid);
newfhdr=fhdr;
newchdr=chdr;

if any(cell_map(:,1) > newfhdr(NCells))
	fclose('all');
	error('Value in source column of cell_map exceeds number of cells in input file.');
end
if any(cell_map(:,2) > newfhdr(NCells))
	fclose('all');
	error('Value in subtractor column of cell_map exceeds number of cells in input file.');
end

newfhdr(NCells)=size(cell_map,1);
newchdr = newchdr(1:size(cell_map,1),:);
newcnames = zeros(size(cell_map,1),80);
for c=1:size(cell_map,1);
	tempstr=['Diff of Cells ' int2str(cell_map(c,1)) ' and ' int2str(cell_map(c,2))];
	newcnames(c,1:length(tempstr)) = tempstr;
end
wt_ave_hdr_v(destfid,newfhdr,newchdr,ename,newcnames,fcom,ftext);

diffdata=zeros(newchdr(1,NSamps),newfhdr(NChan));
for c=1:size(cell_map,1) %c is for cells that live in EGIS
	disp(['Computing difference of cells ' int2str(cell_map(c,1)) ' and ' int2str(cell_map(c,2)) '.'])
	sourcedata = rd_onetr_allch(srcfid, coff(cell_map(c,1)), 1, newfhdr(NChan), chdr(cell_map(c,1), NPoints));
	subtractor = rd_onetr_allch(srcfid, coff(cell_map(c,2)), 1, newfhdr(NChan), chdr(cell_map(c,2), NPoints));
	diffdata = sourcedata - subtractor;
	disp('Writing difference data to output file.');
	fwrite(destfid, diffdata, 'int16');
end	%newcell=1:size(cell_map,1)

disp('Done Running gave_diff.');
fclose('all');
status=1;



