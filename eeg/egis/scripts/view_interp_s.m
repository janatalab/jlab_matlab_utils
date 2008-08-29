function status = view_interp_s(color_map,clims,contour_array,contour_array2,s)
%status = view_interp(color_map,clims,contour_array,contour_array2)
%
%view interpolation binary with header
%
%color_map = what it says
%clims = colormap limits (optional)
%contour_array = first contour array (gets a solid blue contour line) 
%contour_array2 = second contour array (gets a dashed red line). 
fid = get_fid('r');
nrow = fread(fid,1,'float');
ncolumn = fread(fid,1,'float');
plot_res = fread(fid,1,'float');
interp = fread(fid,[ncolumn,nrow],'float');
interp = interp';
if nargin < 5
s = 1;
end;
interp = interp/s;
clims = clims/s;
contour_array = contour_array/s;
contour_array2 = contour_array2/s;
if nargin == 1
clims(2) = 0.9*max(max(interp));
clims(1) = 0.9*min(min(interp));
end;
imagesc(interp,clims), colorbar, colormap(color_map), shading interp
hold on
if (nargin ==2|nargin == 1)
contour(interp)
elseif nargin == 3
contour(interp,contour_array,'r-')
elseif (nargin == 4|nargin == 5)
contour(interp,contour_array,'b-');
contour(interp,contour_array2,'r--');
end;
axis('off')
fclose(fid);
status = 1;



