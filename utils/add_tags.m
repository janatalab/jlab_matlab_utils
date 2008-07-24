function add_tags(filename,tag_prepend,tag_append)
%
% opens a file and adds tag_prepend and tag_append before and after the current
% file contents.
% 
% This is useful for inserting simple html, for example, <table><tr></td> and
% </td></tr></table> tags to the beginning and end of a word or sentence in a
% text file. If you have a series of files to which you need
% to add these tags, you can call this function from a short script that cycles through all of the file names.
%

fid = fopen(filename,'r');
contents = char(fread(fid))';

new_contents = [tag_prepend contents tag_append];
fclose(fid);

fid = fopen(filename,'w');

fwrite(fid,new_contents,'uchar');
fclose(fid);
