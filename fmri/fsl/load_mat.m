function [desmat] = load_mat(fname)
%

unix_str = sprintf('grep -n Matrix %s', fname);
[status,str] = unix(unix_str);

nlines_skip = str2num(strtok(str,':'));

%desmat = textread(fname,'headerlines',nlines_skip+7);
desmat = loadtxt(fname,'skipline',nlines_skip,'convert','force','verbose','off');