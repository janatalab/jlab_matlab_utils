% todd_z.m

LOAD_DATA = 0;
sub_ids = ...
    {'02jun00CW', ...
     '03aug00GG', ...
     '10aug00EH', ...
     '11aug00WS', ...
     '14jun00OMG', ...
     '15jun00JE', ...
     '15jun00SN', ...
     '17jul00RC', ...
     '18jul00CS'}; 

series_map.episodic = [ ...
    3 6; ...
    2 5; ...
    2 6; ...
    3 6; ...
    3 6; ...
    2 5; ...
    3 6; ...
    2 5; ...
    3 6];

series_map.semantic = [ ...
    1 4; ...
    3 6; ...
    3 6; ...
    1 4; ...
    1 4; ...
    3 6; ...
    1 4; ...
    3 6; ...
    1 4];

task_vols = [16:34 50:68 84:102 118:136];
rest_vols = [1:15 35:49 69:83 103:117];

nseries = 6;

if LOAD_DATA 
  load('/home02/janata/data/todd/extracted_data.mat')
end

%
% Create large masks that specify which data to pull from where
%

nvol = size(y,2);
nvol_per_series = nvol/nseries;
nsub = length(sub_ids);

episodic_offsets = (series_map.episodic-1)*nvol_per_series;
semantic_offsets = (series_map.semantic-1)*nvol_per_series;
  
episodic_task_mask = zeros(nvol, nsub)+NaN;
episodic_rest_mask = zeros(nvol, nsub)+NaN;
semantic_task_mask = zeros(nvol, nsub)+NaN;
semantic_rest_mask = zeros(nvol, nsub)+NaN;

for isub = 1:nsub
  idx = task_vols'*ones(1,2) + ones(length(task_vols),1) * episodic_offsets(isub,:);
  episodic_task_mask(idx(:),isub) = 1;

  idx = rest_vols'*ones(1,2) + ones(length(rest_vols),1) * episodic_offsets(isub,:);
  episodic_rest_mask(idx(:),isub) = 1;

  idx = task_vols'*ones(1,2) + ones(length(task_vols),1) * semantic_offsets(isub,:);
  semantic_task_mask(idx(:),isub) = 1;

  idx = rest_vols'*ones(1,2) + ones(length(rest_vols),1) * semantic_offsets(isub,:);
  semantic_rest_mask(idx(:),isub) = 1;
end

%
% Compute z-scores
%

nvox = size(y,1);

for ivox = 1:nvox
  tmp = squeeze(y(ivox,:,:));
  mean_y(ivox,:) = nanmean(tmp.* (episodic_task_mask ...
				  | episodic_rest_mask ...
				  | semantic_task_mask ... 
				  | semantic_rest_mask));
  
  std_y(ivox,:) = nanstd(tmp.* (episodic_task_mask ...
				| episodic_rest_mask ...
				| semantic_task_mask ... 
				| semantic_rest_mask));
  
end

z = zeros(size(y));

for isub = 1:nsub
  z(:,:,isub) = z(:,:,isub) + (y(:,:,isub) - (mean_y(:,isub)* ...
					      ones(1,nvol))) ./ (std_y(:,isub)*ones(1,nvol));
  
end

voxmean = squeeze(mean(y));

epi_task_mean_z = nanmean(voxmean .* episodic_task_mask);
epi_task_std_z = nanstd(voxmean .* episodic_task_mask);

epi_rest_mean_z = nanmean(voxmean .* episodic_rest_mask);
epi_rest_std_z = nanstd(voxmean .* episodic_rest_mask);

sem_task_mean_z = nanmean(voxmean .* semantic_task_mask);
sem_task_std_z = nanstd(voxmean .* semantic_task_mask);

sem_rest_mean_z = nanmean(voxmean .* semantic_rest_mask);
sem_rest_std_z = nanstd(voxmean .* semantic_rest_mask);


datamat = [epi_task_mean_z' epi_rest_mean_z' sem_task_mean_z' ...
	   sem_rest_mean_z'];

task_str = {'epi\_task','epi\_rest','sem\_task','sem\_rest'};