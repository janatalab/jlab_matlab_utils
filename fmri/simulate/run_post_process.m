% run_post_process.m
%
%

% 2003, Petr Janata

des_id = 1;
desmat_dir = fullfile(dataroot, '/test/stats/');

DEMO_STATS = 1;
total_prop_reg  = 0.6; 	% The approximate proportion of variance explained by
                        % the regressors. It will actually be somewhat more
                        % than this

%get_exp_info;
%opt.title = sprintf('%dSil, TR=%1.1f, note_dur=%d ms, cue_stim_soa=%d-%d s, ITI=%1.1f-%1.1f s, total_dur=%2.1f min',tinfo.num_trials(ntt),TR/1000,ms_per_note, min(cue_stim_soa_range)/1000, max(cue_stim_soa_range)/1000, min(inter_trial_interval_range)/1000, max(inter_trial_interval_range)/1000, total_time_min);

%opt.title = strrep(opt.title,'_','\_');
opt.printfig = 0;
opt.append_str = '-append';
%opt.append_str = '';
opt.figstub = fullfile(desmat_dir,sprintf('pb_des%02d', des_id));

% Create a plot that shows the correlations among the regressors
[corrmat] = post_process_corrmat(fullfile(desmat_dir,'SPM_fMRIDesMtx'),opt);

if DEMO_STATS
  
  % Load the design matrix that we created
  load(fullfile(desmat_dir,'SPM_fMRIDesMtx'))
  X = xX.X;
  
  % Remove the constant term that was added by SPM99
  X(:,end) = [];
  
  % Create some data.  
  % Assume that each of the regressors explains an equal amount of variance. 
  
  nreg = size(X,2);
  prop_per_reg = total_prop_reg/nreg
  
  % In creating the fake data, weight some of the regressors of we want to
  weights = ones(1,nreg);
  weights(2:2:end) = 1; % use -1 to make some of the regressors negatively correlated
  
  fakeX = X .* repmat(weights,size(X,1),1);
  
  % Add the regressors together and scale them to between 0 and 1.  This is the
  % part of the fake data than can be explained by our model
  reg_sum = sum(fakeX,2);
  reg_sum = reg_sum-min(reg_sum);
  reg_sum = reg_sum/max(reg_sum);
  
  % Create a noise vector of values between 0 and 1.  This is the part of the
  % data that can't be explained by the model (beyond some chance level)
  rand_vect = rand(size(reg_sum));
  
  % Add the signal and noise together in desired proportion
  fake_data = reg_sum*sqrt(total_prop_reg) +rand_vect*sqrt(1-total_prop_reg);

  % Create a plot that shows the sum of the regressors and the fake data
  figure(3), clf
  subplot(2,1,1)
  plot(reg_sum)
  
  subplot(2,1,2)
  plot(fake_data) 
end

origX = X;

%
% Do some generalized linear model fitting
%

% Fit only a subset of the regressors
use_reg = [1];
X = origX(:,use_reg);
[b_sub, dev, stats_sub] = glmfit(zscore(X), zscore(fake_data), 'normal');
fprintf('Using regressors: %s\n', sprintf('%d ',use_reg));
fprintf('Proportion of variance explained: %s\n', sprintf('%1.2e, ', b_sub(2:end).^2));
fprintf('Significance (p<): %s\n', sprintf('%1.4f, ', stats_sub.p(2:end)))

% Calculate the proportion of variance explained

% Fit the whole model
X = origX;
[b, dev, stats] = glmfit(zscore(X), zscore(fake_data), 'normal');

z = zscore(fake_data);
tss = z'*z;
ess = stats.resid'*stats.resid;

fprintf('\nUsing all regressors\n');
fprintf('Total SS via beta^2: %1.3f\n', sum(b(2:end).^2));
fprintf('Total SS via resid: %1.3f\n', 1-ess/tss);
fprintf('Proportion of variance explained: %s\n', sprintf('%1.2e, ', b(2:end).^2));
fprintf('Significance (p<): %s\n', sprintf('%1.4f, ', stats.p(2:end)))

% Create a plot
bmat = zeros(length(b)-1,2);
pmat = zeros(size(bmat));

bmat(use_reg,1) = b_sub(2:end);
bmat(:,2) = b(2:end);

pmat(use_reg,1) = stats_sub.p(2:end);
pmat(:,2) = stats.p(2:end);

figure(4),clf
subplot(2,1,1)
bar(bmat)
title('Beta values')

subplot(2,1,2)
bar(pmat)
l = line([get(gca,'xlim')],[0.05 0.05]); % draw a red line at a p<0.05 threshold
set(l,'color','r')
title('P values')
