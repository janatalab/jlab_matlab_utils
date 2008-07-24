% plot_intensity.m
%
% Plot timecourse of intensity

dims = size(V);

sum_int_V = sum(reshape(V,prod(dims(1:3)),dims(4)));

pct_diff_max_min = (max(sum_intensity) - min(sum_intensity))/max(sum_intensity)*100;
disp(sprintf('Pct difference between max and min: %1.2f', pct_diff_max_min))

nplot = 5;

subplot(nplot,1,1)
plot(sum_intensity/mean(sum_intensity))

subplot(nplot,1,2)
plot(mean_intensity)

subplot(nplot,1,3)
%intens_rel_mean = (mean_intensity./sum_intensity);
intens_rel_mean = (sum_intensity./mean_intensity);
plot(intens_rel_mean)

subplot(nplot,1,4)
plot(intens_rel_mean/mean(intens_rel_mean))

subplot(nplot,1,5)
vox(1:dims(4)) = V(30,30,20,:);
plot(vox/mean(vox))


  