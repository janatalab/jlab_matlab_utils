function status = view_interp_sc(plot_res,nrows,ncols,divide,scale,color_map,graph_title)
%function status = view_interp_sc(plot_res,nrows,ncols,divide,scale,color_map,graph_title)
%
%Plots and image created by interp_by_samp with the 's'
%header_flag option.

status=0;
fid = get_fid('r');
if fid==-1, error('Don''t pass me no bullshit files, Doran'),end;
interp = fread(fid,[ncols,nrows],'float');
interp = (interp')./divide;
figure;
imagesc(interp, scale);
colorbar;
colormap(color_map);
shading interp;
hold on;
contour(interp);
fclose(fid);
axis('off');
title(graph_title);
status = 1;
