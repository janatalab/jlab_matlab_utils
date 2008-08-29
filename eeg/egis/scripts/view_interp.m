function status = view_interp(color_map,clims,contour_array,contour_array2,title_text,image_labels,export_matrix,nrow,ncolumn,plot_res)
%status = view_interp(color_map,clims,contour_array,contour_array2,image_labelsimage_name,nrow,ncolumn,plot_res)
%
%view interpolation binary with header
%
%color_map = what it says
%clims = colormap limits (optional)
%contour_array = first contour array (gets a solid blue contour line) 
%contour_array2 = second contour array (gets a dashed red line). 
figure
if nargin < 7
[fname,pathname] = uigetfile('*.img','Select image file:');
fid = fopen([pathname fname],'rb');
nrow = fread(fid,1,'float');
ncolumn = fread(fid,1,'float');
plot_res = fread(fid,1,'float');
export_matrix = fread(fid,[ncolumn,nrow],'float');
export_matrix = export_matrix';
fclose(fid);
end
export_matrix(find(export_matrix == -100)) = -0.1*clims(2)*ones(size(find(export_matrix == -100),1),size(find(export_matrix == -100),2))
if nargin == 1
clims(2) = 0.75*max(max(export_matrix));
clims(1) = 1.25*min(min(export_matrix));
end;
if isempty(clims)
clims(2) = 0.6*max(max(export_matrix));
clims(1) = 1.4*min(min(export_matrix));
end;
clims(1) =-0.1*clims(2);
color_map2 = [1 1 1; color_map];
imagesc(export_matrix,clims), colorbar, colormap(color_map2), shading interp
mean_mat = max(max(export_matrix))/4;
steps = 0.1;
max(max(export_matrix))
hold on
if (nargin ==2|nargin == 1)
contour_array2 = [min(min(export_matrix))+step:steps:(mean_mat-steps)];
contour_array = [mean_mat:steps:max(max(export_matrix))];
end;
if isempty(contour_array)
contour_array2 = [(mean_mat-7*steps):steps:(mean_mat-steps)];
contour_array = [mean_mat:steps:max(max(export_matrix))];
end;

if nargin == 3
contour(export_matrix,contour_array,'k-')
else
contour(export_matrix,contour_array,'k-');
contour(export_matrix,contour_array2,'k--');
end;
axis('off')
if nargin >= 5
	title(title_text);
end;
if nargin >= 6
	if size(image_labels,1) == nrow*ncolumn/plot_res/plot_res
		icount = 1;
		nrow = nrow/plot_res;
		ncolumn = ncolumn/plot_res;
		for j = nrow:-1:1
			for i = 1:ncolumn
				text((i-1)*plot_res+fix(plot_res/3),nrow*plot_res-3 -(j-1)*plot_res,image_labels(icount,:))
				icount = icount + 1;
			end
		end
	else
		disp('WARNING:improper number of labels')
		disp('WARNING:improper number of labels')

	end
end

status = 1;












