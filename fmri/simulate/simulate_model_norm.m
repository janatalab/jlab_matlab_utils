% prime_base_model_norm
% 
% This script sets up an event-related analysis.  Although the model for each
% subject is set up from the same prototype, models will differ across
% subjects. One MAJOR DIFFERENCE in models across subjects may be the number
% of conditions modeled.  Since all targets are modeled as a function of
% response accuracy, if subjects made no errors for a particular type of target
% in a particular session, the column modeling that condition will be removed
% from the design matrix.  This means not only that the d.f. will differ across
% subjects but also that the contrasts have to be set up on a subject by
% subject basis.
%


global dataroot spm_path subject_ids use_model_proto sinfo
global SCAN_OFFSET
global REGRESS_MOTION_CORRECTION_PARAMS ADD_LINEAR_TREND_REGRESS ADD_RT_REGRESS

nsub = length(subject_ids);

subj_nsess = length(sinfo(1).cond_order);

disp('Getting subject file info for models')
for isub = 1:nsub
  for isess = 1:subj_nsess
    sessions(isub).images{isess} = spm_get('Files', ... 
	sprintf('%s/%s/epi/run%d/', dataroot, char(subject_ids(isub)), isess), ...
	sprintf('sn%s*.img', char(subject_ids(isub))) ...
	);   
  end  % for isess=
end  % for isub=

prime_long_model_protos;

model_init = model_prototype(use_model_proto);
conditions_init = conditions_prototype(use_model_proto);
 
if REGRESS_MOTION_CORRECTION_PARAMS | ADD_LINEAR_TREND_REGRESS
  regressors_init = regressors_prototype(use_model_proto);
end

ncond_idxs = 0;
nregress_idxs = 0;
  
clear regressors

disp('Getting subject specific model parameters')
for isub = 1:nsub
  model(isub) = model_init;

  model(isub).files = sessions(isub).images;
  
  % Set the SOT values for this subject
  ttinfo_fname = fullfile(dataroot, char(subject_ids(isub)), sprintf('%s_ttinfo.mat', char(subject_ids(isub))));
  sot = prime_SOTs(ttinfo_fname, use_model_proto);
  
  for isess = 1:model(isub).nsess
    ncond_idxs = ncond_idxs + 1;
    sub_conds(isess) = ncond_idxs;		% Keep track of which condition numbers
                                        % correspond to the model for this subject

    conditions(ncond_idxs) = conditions_init;
    conditions(ncond_idxs).onsets = sot{isess};
    
    %
    % Check to make sure that there are events in every condition, and remove
    % conditions for which there are no events
    %
    
    % Compile list of empty conditions
    ncond = length(sot{isess});
    empty_cond = [];
    for icond = 1:ncond
      if isempty(sot{isess}{icond})
	empty_cond(end+1) = icond;
      end
    end
    
    %
    % Delete entries from cond_names, cond_types, onsets, bf_ev, and bf_ep
    %
    
    if empty_cond
      conditions(ncond_idxs).names(empty_cond) = [];
      conditions(ncond_idxs).onsets(empty_cond) = [];
      conditions(ncond_idxs).types(empty_cond) = [];
      conditions(ncond_idxs).bf_ev(empty_cond) = [];
      conditions(ncond_idxs).bf_ep(empty_cond) = [];
      conditions(ncond_idxs).durations{1}(empty_cond) = [];
    end % if empty_cond
    
    sub_nconds(isess) = length(conditions(ncond_idxs).names);

    if REGRESS_MOTION_CORRECTION_PARAMS | ADD_LINEAR_TREND_REGRESS
      nregress_idxs = nregress_idxs + 1;
      regressors(nregress_idxs) = regressors_init;
      nregress_sess = 0;
    end
    
      
    if REGRESS_MOTION_CORRECTION_PARAMS
      % Load regressors
      %fname = spm_get('Files',sprintf('%s/%s/epi/run%d/', dataroot, char(subject_ids(isub)), isess),sprintf('realignment_params_*i%04d.txt', SCAN_OFFSET+1));
      fname = spm_get('Files',sprintf('%s/%s/epi/run%d/', dataroot, char(subject_ids(isub)), isess),sprintf('realignment_params_%s*i%04d.txt', char(subject_ids(isub)),SCAN_OFFSET+1));
      vals = load(fname);
      ncol = size(vals,2);
      regressors(nregress_idxs).values(1:sinfo(isub).nvol(isess)-SCAN_OFFSET,nregress_sess+1:nregress_sess+ncol) = load(fname);
      nregress_sess = nregress_sess + ncol;
    end
    
    if ADD_LINEAR_TREND_REGRESS
      vals = 0:sinfo(isub).nvol(isess)-SCAN_OFFSET-1;
      vals = vals'/max(vals);
      regressors(nregress_idxs).values(1:sinfo(isub).nvol(isess)-SCAN_OFFSET,nregress_sess+1) = vals;
      nregress_sess = nregress_sess + 1;
    end
    
    if ADD_RT_REGRESS
      % load RTs for regressor
      a_info=load(ttinfo_fname); 
      % how many RTs per run?
      nb = length(a_info.rt.raw)/model(isub).nsess;
      a_start = (isess-1)*nb+1; a_stop = isess*nb;
      % for rcr
      a_idx_rcr = find(a_info.good_idx.rcr >= a_start & a_info.good_idx.rcr <= a_stop);
      a_rcr = a_info.good_idx.rcr(a_idx_rcr);
      a_help = a_info.rt.raw(a_rcr);
      a_when = round(sot{isess}{4});
      vals(a_when) = a_help; clear a_help a_when
      % for ucr
      a_idx_ucr = find(a_info.good_idx.ucr >= a_start & a_info.good_idx.ucr <= a_stop);
      a_ucr = a_info.good_idx.ucr(a_idx_ucr);
      a_help = a_info.rt.raw(a_ucr);
      a_when = round(sot{isess}{6});
      vals(a_when) = a_help;clear a_help a_when
      % for rdr
      a_idx_rdr = find(a_info.good_idx.rdr >= a_start & a_info.good_idx.rdr <= a_stop);
      a_rdr = a_info.good_idx.rdr(a_idx_rdr);
      a_help = a_info.rt.raw(a_rdr);
      a_when = round(sot{isess}{5});
      vals(a_when) = a_help;clear a_help a_when
      % for udr
      a_idx_udr = find(a_info.good_idx.udr >= a_start & a_info.good_idx.udr <= a_stop);
      a_udr = a_info.good_idx.udr(a_idx_udr);
      a_help = a_info.rt.raw(a_udr);
      a_when = round(sot{isess}{7});
      vals(a_when) = a_help;clear a_help a_when

      % eliminate the NaN of outliers
      a_nan = isnan(vals); a_wo = find(a_nan > 0); 
      for ai=1:length(a_wo)
	vals(a_wo(ai,1)) = 0;  % eliminate the NaN of outliers
      end 
      
      regressors(nregress_idxs).values(1:sinfo(isub).nvol(isess)-SCAN_OFFSET,nregress_sess+1) = vals;
      nregress_sess = nregress_sess + 1;
    end  

    sub_regress(isess) = nregress_idxs;
    sub_nregress(isess) = length(regressors(nregress_idxs).names);
  end % for isess = 
  
  % Set the condition info for this subject
  model(isub).conditions = sub_conds;
  model(isub).conditions_nb = sub_nconds;
  
  % Set the regressor info for this subject
  model(isub).regressors = sub_regress;  
  model(isub).regressors_nb = sub_nregress;
  
  % Set the correct number of scans for this subject
  model(isub).nscans = sinfo(isub).nvol - SCAN_OFFSET;
  
end  %for isub

