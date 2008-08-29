function [version, NChan, NMeg, NEeg, NReference, NBad_sensor, NTrigger, NResponse, NUtility, NAnalog, Samp_Rate, NData_Epoch, Names, NSamp, trigger, response, header_length] = rd_meg_hdr(meg_fid);
% [version, NChan, NMeg, NEeg, NReference, NBad_Sensor, NTrigger, NResponse,
% NUtility, NAnalog, Samp_Rate, NData_Epoch, Names, NSamp, trigger, response, 
% header_length] = rd_meg_hdr(meg_fid);
%
% reads in header from MEG file.  Note that the file is already open 
% also note that the file is rewinded at the end of this call. 
%
%
%

if meg_fid < 1
	error('you passed rd_hdr_meg a bogus fid')
end;
frewind(meg_fid);
version = fread(meg_fid,1,'integer*2');
if version ~= 1
	error('this is not an meg data file');
end
NChan= fread(meg_fid,1,'integer*2');
NMeg = fread(meg_fid,1,'integer*2');
NEeg = fread(meg_fid,1,'integer*2');
NReference = fread(meg_fid,1,'integer*2');
NBad_sensor = fread(meg_fid,1,'integer*2');
NTrigger = fread(meg_fid,1,'integer*2');
NResponse = fread(meg_fid,1,'integer*2');
NUtility = fread(meg_fid,1,'integer*2');
NAnalog = fread(meg_fid,1,'integer*2');
Samp_Rate = fread(meg_fid,1,'real*4');
NData_Epoch = fread(meg_fid,1,'integer*2');
sizenames = fread(meg_fid,1,'integer*2');
Names = fread(meg_fid,[NChan,sizenames],'char*1');
Names = setstr(Names);
rd_hdr_length = 28 + sizenames*NChan;


NSamp = fread(meg_fid,[1,NData_Epoch],'integer*4');
rd_hdr_length = rd_hdr_length + 4*NData_Epoch;

trigger = fread(meg_fid,[max(NSamp),NTrigger*NData_Epoch],'real*4');

rd_hdr_length = rd_hdr_length + 4*max(NSamp)*NTrigger*NData_Epoch;

response = fread(meg_fid,[max(NSamp), NTrigger*NData_Epoch],'real*4');

rd_hdr_length = rd_hdr_length + 4*max(NSamp)*NResponse*NData_Epoch;

rd_hdr_length = rd_hdr_length + 4;

header_length = fread(meg_fid,1,'integer*4');

if rd_hdr_length ~= header_length
	error('header length dont match up')
end;












