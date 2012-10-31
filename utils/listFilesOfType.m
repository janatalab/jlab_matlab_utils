function flist = listFilesOfType(fpath,ftypes)
% Returns directory listing of all files in diretory specified by fpath of
% one of the types specified in cell-array ftypes

% 27Oct2011 PJ

ntypes = length(ftypes);

flist = {};
for itype = 1:ntypes
	tmplist = dir(fullfile(fpath,sprintf('*.%s',ftypes{itype})));
	
	flist = [flist {tmplist(:).name}];
	
end % for itype

return