function z = fisher_r2z(r);
% z = fisher_r2z(r);
%
% Performs Fisher r to z transformation
%

% July 31, 2009 PJ

z = 0.5*(log(1+r)-log(1-r));

