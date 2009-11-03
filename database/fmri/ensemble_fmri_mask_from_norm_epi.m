function outdata = ensemble_fmri_mask_from_norm_epi(indata,defs)

% mask from normalised epi
% 
% REQUIRES
%   sinfo data - uses to identify subjects and to pass on subject params
%   epi data - takes the first volume of the reference run to create mask -
%       should be normed data, though this script just accepts whatever epi
%       data it gets ... FIXME: should the name of the script be 'from norm
%       epi' or just 'from epi' since it could technically take normed,
%       normed/smoothed, non-normed, etc EPI data??? also, should it just
%       default to placing the data into the defs.fmri.paths.anat_outdir or
%       the epi path if it can't find hires or coplanar data??
%   EITHER hires data OR coplanar data - takes the path to the hires,
%       places the mask on the same path
% 
% FB 2008.08.27

global r

r = init_results_struct;

r.type = 'fmri_anal';  % Identify the type of this reporting instance
r.report_on_fly = 1;

outdata = ensemble_init_data_struct();
outdata.type = 'normed_epi_mask';

% Parse out the input data
for idata = 1:length(indata)
  if isfield(indata{idata},'type')
    switch indata{idata}.type
      case 'sinfo'
        sinfo = indata{idata};
        sinfo = sinfo.data;
        proc_subs = {sinfo(:).id};
        nsub_proc = length(proc_subs);
      case {'epi','realign_epi'}
        epidata = indata{idata};
        epicol = set_var_col_const(epidata.vars);
      case 'hires'
        hires = indata{idata};
        hicol = set_var_col_const(hires.vars);
      case 'coplanar'
        coplanar = indata{idata};
        cocol = set_var_col_const(coplanar.vars);
      case {'paths'}
        pathdata = indata{idata};
        pcol = set_var_col_const(pathdata.vars);
    end
  end
end

% check for required vars
check_vars = {'sinfo','epidata',{'hires','coplanar'}};
check_required_vars;

if (iscell(indata) && ~isempty(indata) && isfield(indata{1},'task') && ...
        ~isempty(strmatch('return_outdir',indata{1}.task))) || ...
        (isstruct(indata) && isfield(indata,'task') && ...
        ~isempty(strmatch('return_outdir',indata.task)))
  if exist('pathdata','var') && ~isempty(pathdata.data{1})
    if length(nsub_proc) == 1
      pfilt = struct();
      pfilt.include.all.subject_id = {proc_subs{1}};
      lpathdata = ensemble_filter(pathdata,pfilt);
      if ~isempty(lpathdata.data{1})
        sfilt = pfilt;
        sfilt.include.all.path_type = {'anat_outdir'};
        spathdata = ensemble_filter(lpathdata,sfilt);
        if length(spathdata.data{1}) == 1
          % one epi outdir, save outdata = epi_outdir
          outdata = spathdata.data{pcol.path}{1};
        else
          sfilt = pfilt;
          sfilt.include.all.path_type = {'sess_outdir'};
          spathdata = ensemble_filter(lpathdata,sfilt);
          if length(spathdata.data{1}) == 1;
            outdata = spathdata.data{pcol.path}{1};
          else
            sfilt = pfilt;
            sfilt.include.all.path_type = {'sub_outdir'};
            spathdata = ensemble_filter(lpathdata,sfilt);
            if length(spathdata.data{1}) == 1;
              outdata = spathdata.data{pcol.path}{1};            
            end
          end
        end
      end
    end
  end
  if ~exist('outdata','var') || ~exist(outdata,'dir'), outdata = ''; end
  return
end

% outdata
outdata.vars = [outdata.vars 'sinfo'];
sinfo_idx = length(outdata.vars);
outdata.data{sinfo_idx} = ensemble_init_data_struct();
outdata.data{sinfo_idx}.type = 'sinfo';
outdata.data{sinfo_idx}.data = sinfo;

outdata.vars = [outdata.vars 'epi_mask'];
emask_idx = length(outdata.vars);
outdata.data{emask_idx} = ensemble_init_data_struct();
outdata.data{emask_idx}.type = 'normed_epi_mask';
outdata.data{emask_idx}.vars = {'subject_id','session','ensemble_id',...
    'path'};
outdata.data{emask_idx}.data{1} = {};
outdata.data{emask_idx}.data{2} = [];
outdata.data{emask_idx}.data{3} = [];
outdata.data{emask_idx}.data{4} = {};

% get flags
try USE_SPM = defs.mask_norm_epi.USE_SPM; catch USE_SPM = 0; end
try USE_FSL = defs.mask_norm_epi.USE_FSL; catch USE_FSL = 0; end

%%%% FIXME: IS THIS CORRECT?
if USE_SPM && ~USE_FSL
  msg = sprintf('SPM not supported yet ...\n');
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
      refFilt.include.all.subject_id = {subid};
      refFilt.include.all.session = isess;
      refFilt.include.all.run = sess.refrun;
      refdata = ensemble_filter(epidata,refFilt);
      
      refvol = 1;
      [rpath,rfname] = fileparts(refdata.data{epicol.path}{refvol});
      srcfname = fullfile(rpath,rfname);

      anatFilt.include.all.subject_id = {subid};
      anatFilt.include.all.session = isess;
      if exist('hires','var')
        anatdata = ensemble_filter(hires,anatFilt);
        [apath] = fileparts(anatdata.data{hicol.path}{1});
      elseif exist('coplanar','var')
        anatdata = ensemble_filter(coplanar,anatFilt);
        [apath] = fileparts(anatdata.data{cocol.path}{1});
      else
        msg = sprintf('no coplanar or hires data found???\n');
        r = update_report(r,msg);
        return
      end
      
      % Use bet2 to do the mask creation
      maskstub = fullfile(apath,sprintf('sw%s_%s', subid, sess.id));
      unix_str = sprintf('bet2 %s %s -m',srcfname, maskstub);
      msg = sprintf('%s\n', unix_str);
      r = update_report(r,msg);
      unix(unix_str);
      
      % gunzip the mask that was just created, or SPM will croak      
      gunstr = sprintf('%s_mask.nii',maskstub);
      gstr = sprintf('%s.gz',gunstr);
      if exist(gstr,'file')
        if exist(gunstr,'file')
          delete(gunstr);
        end
        unix_str = sprintf('gunzip %s',gstr);
        msg = sprintf('%s\n',unix_str);
        r = update_report(r,msg);
        unix(unix_str);
      elseif ~exist(gunstr,'file')
        error('mask file (%s) not created!?\n',gunstr);
      end
      
      outdata.data{emask_idx}.data{1} = [outdata.data{emask_idx}.data{1} ...
          subid];
      outdata.data{emask_idx}.data{2} = [outdata.data{emask_idx}.data{2} ...
          isess];
      outdata.data{emask_idx}.data{3} = [outdata.data{emask_idx}.data{3} ...
          sess.ensemble_id];
      outdata.data{emask_idx}.data{4} = [outdata.data{emask_idx}.data{4} ...
          sprintf('%s_mask.nii',maskstub)];
    end % if USE_FSL
    
  end % for isess
end % for isub=
