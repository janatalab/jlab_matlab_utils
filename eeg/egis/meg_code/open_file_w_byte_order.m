function [fid] = open_file_w_byte_order(filename,version_parameter);
% [fid] = open_file_w_byte_order(filename,version_parameter);
%
%
fid = fopen(filename,'rb','n');
version = fread(fid,1,'int16');

if version ~= version_parameter
	fclose(fid);
	fid = fopen(filename,'rb','b');
	version = fread(fid,1,'int16');
end;

if version ~= version_parameter
	fclose(fid);
	fid = fopen(filename,'rb','l');
	version = fread(fid,1,'int16');
end;

if version ~= version_parameter
	error('tried both byte orders. this is not the right type of file');
end;

frewind(fid);