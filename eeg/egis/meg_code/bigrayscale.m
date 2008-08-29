function cmap = bigrayscale(N)
% Produce antisymmetric grayscale cmap for displaying signed data in
% conjunction with contours.  Has 2N-1 shades, ranging from white to white,
% passing through black in the middle.

if ~exist('N'), N=64; end
a=gray(N);
b=flipud(a);
cmap = [b(1:(N-1),:); a];