% counter_sim.m
%
% Illustrates that there is no theoretical difference between counter.m and
% randperm.  Randperm is faster
%
%

% 11/24/03  Petr Janata

ntrials = 20;
ncond = 4;
total_trials = ncond*ntrials;

% Create a list of 4 conditions, 20 trials each
fprintf('Creating design matrix with randperm ...\n');
X = reshape(ones(ntrials,1)*(1:ncond),total_trials,1);

% Create a dummy coded matrix with conditions in columns.  Randomize assignment
% of condition to each trial using randperm
X = full(sparse(1:total_trials,X(randperm(total_trials)),1));

% Eliminate the 4th condition
X(:,ncond) = [];

% Create a similar matrix using counter.m for the randomization
fprintf('Creating desing matrix with counter ...\n');
[c] = counter(ncond, ntrials*ncond);

Y = full(sparse(1:total_trials,c,1));
Y(:,ncond) = [];

% Now dump the correlation matrix for each design matrix

fprintf('Correlation matrix generated using randperm:\n')
corrcoef(X)

fprintf('\nCorrelation matrix generated using counter:\n')
corrcoef(Y)