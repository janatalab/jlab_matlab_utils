function z = ppc_cmplx_mult(x,y)
% function z = ppc_cmplx_mult(x,y)
%
% multiplies two complex valued vectors by breaking up
% the multiplication into multiple steps to avoid rounding
% errors that will return erroneous results if a number is multiplied
% by its complex conjugate on the PowerMac
%
% Use this function if you want a complex vector multiplied by its conjugate
% to return a real-valued vector.
%
% Tech-support at Mathworks suggested this fix.

% Written by Petr Janata, Jan. 17, 1996.

z = real(x).*real(y) + real(x).*imag(y) + real(y).*imag(x) - imag(x).*imag(y);

return
