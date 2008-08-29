function [avgimg] = average_interp(interp_data,plot_res,samps2average)
% average_interp.m
%
% interp_data can either be a plot_res X (plot_resXnsamps) image, or a plot_res
% X plot_res X nsamps image
%

% 06/13/00  PJ  Modified to accomodate interp_data of various sizes

dims = size(interp_data);
ndim = length(dims);

if ndim == 3
  if (dims(1) ~= dims(2)) & (dims(1) ~= plot_res)
    error('Do not understand dimensionality: %s', sprintf('%d ', dims))
  end
end

navg = length(samps2average);

avgimg = 0;

if ndim == 2
  for isamp = 1:navg
    avgimg = avgimg + ...
	interp_data(:,(samps2average(isamp)-1)*plot_res+1:samps2average(isamp)*plot_res);
  end
  avgimg = avgimg/navg;
else
  avgimg = mean(interp_data(:,:,samps2average),3);
end

return
