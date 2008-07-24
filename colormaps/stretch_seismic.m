function ns = stretch_seismic(map,stretch_factor)
%
% Script not finished

% Find index at which green goes to 1

gidx = find(map(:,2) > 0 & map(:,2) < 1);

grbrk(1) = gidx(find(diff(gidx))>1);
grbrk(2) = gidx(find(diff(diff(gidx))<0)+1);

ns = map;


