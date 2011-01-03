function [names,vals] = fmri_regress_cardiac_liston(pinfo,minfo,sess)

% generates cardiac regressors based on phase, one reg per slice
% 
%   [names,vals] = fmri_regress_cardiac_liston(pinfo,minfo,sess)
% 
% TAP method, Liston et al (2006) NeuroImage
% 
% the TAP method published by Liston consisted of a single phase
% calculation for each slice, for each volume, and it includes a
% separate regressor for each slice in a volume-based model.
%.Dagli e.a. (1999), Glover e.a. (2000), and others have used a
% 6-dimension sin and cos transformation of a single phase
% calculation for each slice/volume, in a slice-based model.
% Here we give the option of expanding the Liston volume-based model
% to include 3 sin and 3 cos transformations from Dagli and Glover's
% slice-based models within a volume-based model
% 
% called by fmri_generate_regress
% 
% REQUIRES
% 
% RETURNS
%   names = cell array containing one regressor name: 'linear_run%d'
%   vals = -(log(t/T)/log(t(1)/T)) + 1
% 
% FB 2009.11.05

% init vars
TR = pinfo.scanner.TR;
nslice_per_vol = pinfo.scanner.orig_nslices;
nvol = pinfo.scanner.actual_nvol;

% init output vars
names = cell(1,nslice_per_vol*6);
vals = zeros(nvol,nslice_per_vol*6);

% get cardiac data
card_epochs = pinfo.card_epochs;
ccol = set_var_col_const(card_epochs.vars);
pidxs = card_epochs.data{ccol.peakidxs}{1};
srate = card_epochs.data{ccol.srate}(1);

if iscell(pidxs), pidxs = pidxs{1}; end
sidxs = pidxs./srate; % peak indices in seconds

% get slice onsets
try slice_order_type = pinfo.protocol.slice_acquisition_order;
catch slice_order_type = 'interleaved'; end
slice_dur = TR/nslice_per_vol;
slice_times = 0:slice_dur:(TR-slice_dur);

sonset_mtx = zeros(nvol,nslice_per_vol);

switch slice_order_type
    case 'interleaved'
      slice_order = [1:2:nslice_per_vol 2:2:nslice_per_vol];
      [sorted,unsorting] = sort(slice_order);
      sorted_slice_times = slice_times(unsorting);
    case 'sequential'
      sorted_slice_times = slice_times;
    otherwise
      error('unknown slice order type: %s',slice_order_type);
end

% generate volume X slice onset time matrix
for ivol = 1:nvol
  sonset_mtx(ivol,:) = sorted_slice_times + (ivol-1)*TR;
end

% find phase of cardiac signal at each slice
for ions=1:numel(sonset_mtx)

  c1i = find(sidxs <= sonset_mtx(ions),1,'last');
  c2i = find(sidxs > sonset_mtx(ions),1,'first');

  if isempty(c1i), c1 = 0; else c1 = sidxs(c1i); end
  if isempty(c2i)
    % estimate last cardiac pulse in the cycle
    c2 = sidxs(end) + mean(diff(sidxs(end-9:end)));
  else
    c2 = sidxs(c2i);
  end
        
  % As described in Dagli et al. 1999, NeuroImage 9, 407-415
  sonset_mtx(ions) = (sonset_mtx(ions) - c1)/(c2-c1);
end

% sin transforms
for is=1:3
  sonset_mtx(:,:,end+1) = sin(2*pi*sonset_mtx(:,:,1)*is);
end

% cos transforms
for is=1:3
  sonset_mtx(:,:,end+1) = cos(2*pi*sonset_mtx(:,:,1)*is);
end

sonset_mtx(:,:,1) = [];

stub6 = {'sin1','sin2','sin3','cos1','cos2','cos3'};

for is=1:numel(sonset_mtx(1,:,:))
  snum = mod(is-1,nslice_per_vol)+1;
  rnum = fix((is-1)/nslice_per_vol)+1;

  names{is} = sprintf('cardiac_%s_slice%1.0f',stub6{rnum},snum);

  vals(:,is) = squeeze(sonset_mtx(:,is));
end % for is=1:prod(size(tap_mtx
