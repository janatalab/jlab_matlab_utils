function cos_array = ang_dist(sources,electrodes)
% cos_array = ang_dist(sources,electrodes)
% sources = theta,phi 
% electrodes =theta,phi
cos_array = zeros(size(sources,1),size(electrodes,1));

d2_t_s = sources(:,1)*ones(1,size(electrodes,1));
d2_p_s = sources(:,2)*ones(1,size(electrodes,1));

d2_t_e = ones(size(sources,1),1)*(electrodes(:,1)');
d2_p_e = ones(size(sources,1),1)*(electrodes(:,2)');

cos_array = cos(d2_t_e).*cos(d2_t_s)+sin(d2_t_e).*sin(d2_t_s).*cos(d2_p_e - d2_p_s);
