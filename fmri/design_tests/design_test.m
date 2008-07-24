% design_test.m

global rootpath

rootpath = '/data1/design_tests/';

%
%  Create data vector
%

means = [5 2 8 1];
ncond = length(means);
nruns = 2;

npts_per_cond = 20;

npts = ncond*nruns*npts_per_cond;

data = repmat(ones(npts_per_cond,1)*means,[1 1 nruns]);
data = data(:);

noise_scaling = 4;
noise = rand(size(data));
noise = noise - mean(noise);
noise = noise*noise_scaling;

noisy = noise+data;

plot(noisy)

%
%  Create images
%

spm_bch('design_batch');

load SPM_fMRIDesMtx

