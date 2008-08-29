function [status]=wt_csdm_hdr_v(fid,infhdr,inchdr,ename,cnames,fcom,ftext);
%[status]=wt_csdm_hdr_v(fid,fhdr,chdr,ename,cnames,fcom,ftext)
%
%
% 1/10/96 PJ Changing calculation of LData to reflect fact that 4-byte floats
%         are now being written
%
% 1/9/96  PJ & RS Length of Header is being written into old BytOrd field
%         This gives it int32 representation to handle overly large header.
%
%This function creates a new EGIS average file header in the file with id fid.
%The structures which are passed to this routine are the same as
%are returned from function rd_ses_hdr_v (c.f.) and which are
%passed to wt_ses_hdr_v.  
%This routine reformats the file header and recalculates its size.
%The total number of elements written is returned.
%Note that this number is NOT the header size in bytes.
%
%NOTE: 	The NObs (NTrials) and NAvg (number of trials in average) fields in 
%		the Average Cell Specifics must be set prior to calling this routine
%		The same if true for the LastDone field of fhdr (Average File Description)
%
%See also:rd_ses_hdr_v.m

%  Modification history:
%	6/6/95	PJ	Created the file based on wt_ses_hdr_v
%

size_data_value = 4;  % 4=float

csdm_hdr_offsets_v;

fhdr = infhdr;
chdr = inchdr;

%
% Perform necessary transformations on the header
%

%% Modify Average File Description
	fhdr(HdrVer) = -1;
	
%% Modify Average Information
	fhdr(5:10) = zeros(1,6);
	fhdr(ScaleBins) = 500;
	fhdr(ScaleCal) = 50;
	fhdr(BaseDur) = 0;
	fhdr(15:17) = zeros(1,3);

%% Modify Average File Data Description
	if (fhdr(NChan) == 64 | fhdr(NChan) == 128)
		fhdr(NChan) = fhdr(NChan) + 1;	% add a channel now that the average reference has been computed
	end
	
%% Modify Average Cell Headers
if any(chdr(:, LSpec) < 28)
  error(['Length of specs in cell(s) [' find(chdr(:,LSpec) < 28) '] is less than required minimum']);
end

if infhdr(HdrVer) == 1	% session file
	ses_hdr_offsets_v;
	chdr(:,RunDateMo_Ave) = ones(fhdr(NCells),1).*infhdr(RunDateMo);
	chdr(:,RunDateDay_Ave) = ones(fhdr(NCells),1).*infhdr(RunDateDay);
	chdr(:,RunDateYr_Ave) = ones(fhdr(NCells),1).*infhdr(RunDateYr);
	chdr(:,RunTimeHr_Ave) = ones(fhdr(NCells),1).*infhdr(RunTimeHr);
	chdr(:,RunTimeMin_Ave) = ones(fhdr(NCells),1).*infhdr(RunTimeMin);
	chdr(:,RunTimeSec_Ave) = ones(fhdr(NCells),1).*infhdr(RunTimeSec);
	chdr(:,SubjID_Ave) = ones(fhdr(NCells),1).*infhdr(SubjID);
	chdr(:,Handed_Ave) = ones(fhdr(NCells),1).*infhdr(Handed);
	chdr(:,Sex_Ave) = ones(fhdr(NCells),1).*infhdr(Sex);
	chdr(:,Age_Ave) = ones(fhdr(NCells),1).*infhdr(Age);
	chdr(:,ExperID_Ave) = ones(fhdr(NCells),1).*infhdr(ExperID);
	chdr(:,EdVer_Ave) = ones(fhdr(NCells),1).*infhdr(EdVer);
	chdr(:,CalFlag_Ave) = ones(fhdr(NCells),1);
end

%% Recompute header length

%disp('Computing length variables in header of average EGIS file...')
L = ((chdr(:,LSpec).*chdr(:,NObs))'+90);
fhdr(LCellHdr:(LCellHdr+(fhdr(NCells)-1))) = L';
nlcellhdr = fhdr(NCells)*2;
lhdr = 130+nlcellhdr+ fhdr(NChan) * 4 + sum(L)+fhdr(LComment)+fhdr(LText);
fhdr(LPad) = (512-rem(lhdr,512));
lhdr = lhdr+fhdr(LPad);
fhdr(LHeader_CSDM) = lhdr;
fhdr(LData) = sum((chdr(:,NPoints).*chdr(:,NObs))*fhdr(NChan))*size_data_value;

status=fseek(fid,0,'bof');
status=fwrite(fid,fhdr(LHeader_CSDM),'int32');
message=ferror(fid);
if message ~= '',error(message),end

status=fwrite(fid,fhdr(HdrVer),'int16');
message=ferror(fid);
if message ~= '', error(message),end

status=fwrite(fid,fhdr(LHeader),'int16');
message=ferror(fid);
if message ~= '',error(message),end

status=fwrite(fid,fhdr(LData),'int32');
message=ferror(fid);
if message ~= '',error(message),end

status=fwrite(fid,abs(ename),'char');
message=ferror(fid);
if message ~= '',error(message),end

status=fwrite(fid,zeros(6,1),'int16');  % empty fields -- run date and time in ses files
message=ferror(fid);
if message ~= '',error(message),end

status=fwrite(fid,fhdr(LastDone:BaseDur),'int16');
message=ferror(fid);
if message ~= '',error(message),end

status=fwrite(fid,zeros(3,1),'int16');
message=ferror(fid);
if message ~= '',error(message),end

status=fwrite(fid,fhdr(NCells:BrdGain),'int16');
message=ferror(fid);
if message ~= '',error(message),end

status=fwrite(fid,fhdr(LCellHdr:(LCellHdr+(fhdr(NCells)-1))),'int16');
message=ferror(fid);
if message ~= '',error(message),end

status=fwrite(fid,zeros(fhdr(NChan)*2, 1),'int16');  % fill with zeros
message=ferror(fid);
if message ~= '',error(message),end

for loop=1:fhdr(NCells)
	status=fwrite(fid,chdr(loop,CellID),'int16');
	message=ferror(fid);
	if message ~= '',error(message),end
	status=fwrite(fid,abs(cnames(loop,:)),'char');
	message=ferror(fid);
	if message ~= '',error(message),end
	lastrow=5+((chdr(loop,LSpec)/2)*chdr(loop,NObs));
	status=fwrite(fid,chdr(loop,NObs:lastrow),'int16');
	message=ferror(fid);
	if message ~= '',error(message),end
end

status=fwrite(fid,fcom,'char');
message=ferror(fid);
if message ~= '',error(message),end

status=fwrite(fid,ftext,'char');
message=ferror(fid);
if message ~= '',error(message),end

status=fwrite(fid,zeros(fhdr(LPad),1),'char');
message=ferror(fid);
if message ~= '',error(message),end

status = 0;
disp(['Successfully wrote the average file header. ' '(' num2str(fhdr(LHeader_CSDM)) ' bytes long)']);
