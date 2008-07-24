function cmap_idx = data2cmap(data,map_size,clim)
% cmap_idx = data2cmap(data,map_size,prop_abs);
%
% Converts data into colormap indices
%
% map_size is the number of rows in the colormap (default: 64)
% clim gives the endpoints of the scale (default: [min(data(:)) max(data(:))])

% 12/1/04 PJ

try map_size(1); catch map_size = 64; end;
try clim(1); catch clim = []; end;

% Determine the colormap scaling
if isempty(clim)
  cscale = linspace(min(data(:)), max(data(:)), map_size);
else
  cscale = linspace(clim(1),clim(2),map_size);
end

[count,bin_idxs] = histc(data(:),cscale);

% Deal with out of range values
bin_idxs(data < min(cscale)) = 1;
bin_idxs(data > max(cscale)) = map_size;

cmap_idx = reshape(bin_idxs,size(data));

return