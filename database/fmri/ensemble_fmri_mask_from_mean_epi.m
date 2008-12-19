function outdata = ensemble_fmri_mask_from_mean_epi(indata,defs)

% creates mask from mean epi
% 
% REQUIRES
%   sinfo
%   epi data
%   mean_epi data
% 
% FB 2008.08.27

global defs r

r = init_results_struct;

r.type = 'fmri_anal';  % Identify the type of this reporting instance
r.report_on_fly = 1;

outdata = ensemble_init_data_struct();
outdata.type = 'mean_epi_mask';

% Parse out the input data
for idata = 1:length(indata)
  switch indata(idata).type
    case 'sinfo'
      sinfo = indata(idata);
    %%%% FIXME: get mean EPI data ... where would that come from???
    case {'mean_epi'}
      medata = indata(data);
      mecol = set_var_col_const(medata.vars);
  end
end

% check for required vars
good = true;
check_vars = {'sinfo','medata'};
for icv=1:length(check_vars)
  if ischar(check_vars{icv})
    if ~exist(check_vars{icv},'var')
      msg = sprintf('\nCOULD NOT FIND VALID SINFO DATA STRUCT\n');
      r = update_report(r,msg);
      good = false;
    end
  elseif iscell(check_vars{icv}) && good
    % check that at least one of the vars exists, otherwise ~good
    conditional = check_vars{icv};
    good=false;
    for iccv=1:length(conditional)
      if exist(conditional{iccv},'var')
        good = true;
      end
    end
  end
end

if ~good, return, end

sinfo = sinfo.data;
proc_subs = {sinfo(:).id};

% outdata
outdata.vars = [outdata.vars 'sinfo'];
sinfo_idx = length(outdata.vars);
outdata.data{sinfo_idx} = ensemble_init_data_struct();
outdata.data{sinfo_idx}.type = 'sinfo';
outdata.data{sinfo_idx}.data = sinfo;

outdata.vars = [outdata.vars 'epi_mask'];
emask_idx = length(outdata.vars);
outdata.data{emask_idx} = ensemble_init_data_struct();
outdata.data{emask_idx}.type = 'mean_epi_mask';
outdata.data{emask_idx}.vars = {'subject_id','session','ensemble_id',...
    'path'};
outdata.data{emask_idx}.data{1} = {};
outdata.data{emask_idx}.data{2} = [];
outdata.data{emask_idx}.data{3} = [];
outdata.data{emask_idx}.data{4} = {};

% get flags
try USE_SPM = defs.mask_mean_epi.USE_SPM; catch USE_SPM = 0; end
try USE_FSL = defs.mask_mean_epi.USE_FSL; catch USE_FSL = 0; end

%%%% FIXME: IS THIS CORRECT?
if USE_SPM && ~USE_FSL
  msg = sprintf('FSL not supported yet ...\n');
  r = update_report(r,msg);
  return
elseif ~USE_FSL && ~USE_SPM
  msg = sprintf(['\t\tyou must specify either SPM or FSL to carry out '...
      'the analyses (but for this analysis, only FSL is currently supported)\n']);
  r = update_report(r,msg);
  return
end

%
% START OF THE SUBJECT LOOP
%

nsub_proc = length(proc_subs);

for isub=1:nsub_proc
  subid = sinfo(isub).id;
  msg = sprintf('\t\tPROCESSING SUBJECT (%d/%d): %s\n', isub, nsub_proc,subid);
  r = update_report(r,msg);

  % Determine number of sessions for this subject
  nsess = length(sinfo(isub).sessinfo);
  
  %
  % START OF THE SESSION LOOP
  %
  
  for isess = 1:nsess
    sess = sinfo(isub).sessinfo(isess);
    
    if ~sess.use_session
      msg = sprintf('\t\t\tSkipping session %d\n', isess);
      r = update_report(r,msg);
      continue
    end

    if USE_FSL
      refFilt.include.all.subject_id = subid;
      refdata = ensemble_filter(medata,refFilt);

      [rpath,rfname] = fileparts(refdata.data{mecol.path}{1});
      srcfname = fullfile(rpath,rfname);

      % Use bet2 to do the mask creation
      maskstub = fullfile(rpath,sprintf('%s', subid));
      unix_str = sprintf('bet2 %s %s -m',srcfname, maskstub);
      msg = sprintf('%s\n', unix_str);
      r = update_report(r,msg);
      unix(unix_str);

      outdata.data{emask_idx}.data{1} = [outdata.data{emask_idx}.data{1} ...
          subid];
      outdata.data{emask_idx}.data{2} = [outdata.data{emask_idx}.data{2} ...
          isess];
      outdata.data{emask_idx}.data{3} = [outdata.data{emask_idx}.data{3} ...
          sess.ensemble_id];
      outdata.data{emask_idx}.data{4} = [outdata.data{emask_idx}.data{4} ...
          sprintf('%s_mean.nii',maskstub)];
    end % if USE_FSL

  end % for isess=
end % for isub=
