function [status]=wt_ave_hdr_v(fid,infhdr,inchdr,ename,cnames,fcom,ftext);
%[status]=wt_ave_hdr_v(fid,fhdr,chdr,ename,cnames,fcom,ftext)
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
%	1/12/96 PJ	BytOrd value is set to 16909060 (integer*4) so that
%			header modified by wt_csdm_hdr_v is reset to default
%			value
%
%	6/6/95	PJ	Created the file based on wt_ses_hdr_v
%
% 04/30/97 PJ  Fixed machine architecture compatability

if nargin < 7
  error('Too few arguments to wt_ave_hdr_v')
end

ave_hdr_offsets_v;
fhdr = infhdr;
chdr = inchdr;

%
% Perform necessary transformations on the header
%

%% Put in the proper BytOrd value -- this is changed by wt_csdm_hdr_v, 
%% so it has to be changed back when writing files derived from csdm files

	fhdr(BytOrd) = 16909060;

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
	ses_hdr_offsets;
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

% Check to make sure that cell name matrix has 80 columns
if size(cnames,2) ~= 80
  disp('Resizing cell name array')
  if size(cnames,2) < 80
    cnames = [cnames zeros(size(cnames,1),80-size(cnames,2))];
  else
    cnames = cnames(:,1:80);
  end
end

%% Recompute header length
fhdr(LComment) = length(fcom);
fhdr(LText) = length(ftext);


%disp('Computing length variables in header of average EGIS file...')
L = ((chdr(:,LSpec).*chdr(:,NObs))'+90);
fhdr(LCellHdr:(LCellHdr+(fhdr(NCells)-1))) = L';
nlcellhdr = fhdr(NCells)*2;
lhdr = 130+nlcellhdr+ fhdr(NChan) * 4 + sum(L)+fhdr(LComment)+fhdr(LText);
fhdr(LPad) = (512-rem(lhdr,512));
lhdr = lhdr+fhdr(LPad);
fhdr(LHeader) = lhdr;
fhdr(LData) = sum((chdr(:,NPoints).*chdr(:,NObs))*fhdr(NChan))*2;

status=fseek(fid,0,'bof');
status=fwrite(fid,fhdr(BytOrd),'integer*4');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,fhdr(HdrVer),'integer*2');
message=ferror(fid);
if ~isempty(message), error(message),end

status=fwrite(fid,fhdr(LHeader),'integer*2');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,fhdr(LData),'integer*4');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,abs(ename),'char');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,zeros(6,1),'integer*2');  % empty fields -- run date and time in ses files
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,fhdr(LastDone:BaseDur),'integer*2');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,zeros(3,1),'integer*2');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,fhdr(NCells:BrdGain),'integer*2');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,fhdr(LCellHdr:(LCellHdr+(fhdr(NCells)-1))),'integer*2');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,zeros(fhdr(NChan)*2, 1),'integer*2');  % fill with zeros
message=ferror(fid);
if ~isempty(message),error(message),end

for loop=1:fhdr(NCells)
	status=fwrite(fid,chdr(loop,CellID),'integer*2');
	message=ferror(fid);
	if ~isempty(message),error(message),end
	status=fwrite(fid,abs(cnames(loop,:)),'char');
	message=ferror(fid);
	if ~isempty(message),error(message),end
	lastrow=5+((chdr(loop,LSpec)/2)*chdr(loop,NObs));
	status=fwrite(fid,chdr(loop,NObs:lastrow),'integer*2');
	message=ferror(fid);
	if ~isempty(message),error(message),end
end

status=fwrite(fid,fcom,'char');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,ftext,'char');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,zeros(fhdr(LPad),1),'char');
message=ferror(fid);
if ~isempty(message),error(message),end

status = 0;
disp(['Successfully wrote the average file header. ' '(' num2str(fhdr(LHeader)) ' bytes long)']);
