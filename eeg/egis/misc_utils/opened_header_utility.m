% rd_ses_hdr.m
%
% Reads the header information in EGIS session (.ses) files

% Modification history:
%
%  2/15/95 PJ  Created this m-file from the original egis.m file
%

[fid, message]=fopen(fname,'rb','b');
if fid == -1
  message
end

offset = 6;
status = fseek(fid, offset, 'bof');
% read LHeader
LHeader = fread(fid, 1, 'int16');
LData = fread(fid, 1, 'int32');
ExpName = fread(fid, 80, 'char');

offset = 12;
status = fseek(fid, offset, 'cof');

SubjID = fread(fid, 1, 'int16');
Handed = fread(fid, 1, 'int16');
Sex = fread(fid, 1, 'int16');
Age = fread(fid, 1, 'int16');
ExperID = fread(fid, 1, 'int16');

offset = 4;
status = fseek(fid, offset, 'cof');

NCells = fread(fid, 1, 'int16');
NChan = fread(fid, 1, 'int16');

offset = 8;
status = fseek(fid, offset, 'cof');

LCellHeader = fread(fid, NCells, 'int16');
CalZero = fread(fid, NChan, 'int16');
CalGain = fread(fid, NChan, 'int16');

%Read in the Cell Header data for all the cells
Cell_ID = zeros(1,NCells);
CellName = zeros(NCells, 80);
Ntrials = zeros(1, NCells);
Npoints = zeros(1, NCells);
SampRate = zeros(1, NCells);
LSpec = zeros(1, NCells);
TrialSpecs = zeros(Ntrials(1)*64, NCells);

DatumSize = 2;
Cell_DataOffset = zeros(1,NCells);	% Calculate beginning of each
					% cell's data

cum_cell_sizes = 0;

for c = 1:NCells
  Cell_ID(c) = fread(fid, 1, 'int16');
  offset = ftell(fid) + 80;
  fscanf(fid, '%s', 80);
  fseek(fid, offset, 'bof');
  Ntrials(c) = fread(fid, 1, 'int16');
  Npoints(c) = fread(fid, 1, 'int16');
  SampRate(c) = fread(fid, 1, 'int16');
  LSpec(c) = fread(fid, 1, 'int16');
  LSpec(c) = LSpec(c)/2;

  TrialSpecs(1:LSpec(c).*Ntrials(c), c) = fread(fid, [LSpec(c).*Ntrials(c),1], 'int16');

  Cell_DataOffset(c) = LHeader + cum_cell_sizes;
  cum_cell_sizes = cum_cell_sizes + Ntrials(c)*Npoints(c)*NChan*DatumSize;
end

fseek(fid, LHeader, 'bof');











