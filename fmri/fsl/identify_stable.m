function sc = identify_stable(cdata)
% sc = identify_stable(corr_data);
%
% Identifies paths of correlations through the correlation matrices in cdata
%

% 10/19/05 Petr Janata - extracted from earlier script

nsess = size(cdata,1);

for isess = 1:nsess
  tmp = sparse(1);
  for jsess = 1:nsess
    sess1 = mod(isess+jsess-2,nsess)+1;
    sess2 = mod(sess1,nsess)+1;
    tmp = tmp*cdata{sess1,sess2};
  end
  sc{isess} = diag(tmp);
end % for isess=