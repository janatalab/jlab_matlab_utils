function corr_data = corr_ics(flist, ccrit);
% corr_data = corr_ics(flist, ccrit);
%
% Correlates Independent Components (or any image volumes for that matter) in
% the list of 4D image files that are given in the cell string variable flist.
% ccrit is the correlation coefficient criterion.
%
% If flist consists of data for nsess sessions,
% then corr_data is an nsess X nsess cell array in which each element is a
% sparse matrix of size invol X jnvol, where invol and jnvol are the sizes of
% the 4th dimensions (ICs or time) of the ith and jth sessions in flist.  Only
% those correlations exceeding ccrit are stored in the sparse matrices.
%

% 10/19/05 Petr Janata - extracted from an earlier script for doing this

nsess = length(flist);

for isess = 1:nsess
  isess_fname = flist{isess};

  % Get the size of the 4th dimension
  [dummy, tmpstr] = unix(sprintf('avwnvols %s', isess_fname));
  invol = str2num(tmpstr);
  
  for jsess = isess:nsess
    jsess_fname = flist{jsess};
    
    % Get the size of the 4th dimension    
    [dummy, tmpstr] = unix(sprintf('avwnvols %s', jsess_fname));
    jnvol = str2num(tmpstr);
    
    unix_str = sprintf('avwcc %s %s %1.2f', isess_fname, jsess_fname, ccrit);
    fprintf('%s\n', unix_str);
    [dummy, tmpstr] = unix(unix_str);
    [iidx, jidx, ccval] = strread(tmpstr(1:end-1),'%d%d%f','delimiter',' ');
    corr_data{isess,jsess} = sparse(iidx,jidx,ccval,invol,jnvol);
    corr_data{jsess,isess} = sparse(jidx,iidx,ccval,jnvol,invol);
  end % for jsess
end % for isess

return