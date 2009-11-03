function outdata = ensemble_fmri_apply_fsl_xfm(indata,defs)

% applies FSL registration parameters to the given EPI files
% 
%   outdata = tarp_fmri_calc_physio_var_subject(indata,defs)
% 
% Given a set of 4D volumes and subject info, this script tries to find,
% and then apply, FSL registration transformations to the 4D volumes. If a
% sequence is found, such as epi2coplanar, coplanar2hires, hires2template,
% this script uses these transformation files to compute epi2template. It
% uses the coregistration formula found in defs.
% 
% REQUIRES
%   indata
%       sinfo - ensemble data struct containing sinfo for all subjects to
%           be processed
%       epidata/epi/realign_epi
%       resid - ensemble data struct containing paths to residuals for each
%           model in feat_dirs
%       rsqrpaths
%   defs
%       defs.fmri.
% 
% RETURNS
% 
% FB 2009.10.30

% % % INIT VARIABLES
outdata = '';
dstr = datestr(now(),30);

% return outdir?
if (iscell(indata) && ~isempty(indata) && isfield(indata{1},'task') && ...
        ~isempty(strmatch('return_outdir',indata{1}.task))) || ...
        (isstruct(indata) && isfield(indata,'task') && ...
        ~isempty(strmatch('return_outdir',indata.task)))
  return
end

display('initializing parameters and paths');

outdata = ensemble_init_data_struct();
outdata.type = 'apply_fsl_xfm';

r = init_results_struct;

r.type = 'apply_fsl_xfm';  % Identify the type of this reporting instance
r.report_on_fly = 1;

% Parse out the input data
for idata = 1:length(indata)
  if isfield(indata{idata},'type')
      switch indata{idata}.type
        case 'sinfo'
          sinfo = indata{idata}.data;
        case {'epi','epidata','realign_epi'}
          epidata = indata{idata};
          epicol = set_var_col_const(epidata.vars);
        case {'resid'}
          residata = indata{idata};
          rescol = set_var_col_const(residata.vars);
        case 'rsqrpaths'
          rsqrpaths = indata{idata};
          rscol = set_var_col_const(rsqrpaths.vars);
        case 'maskpaths'
          maskpaths = indata{idata};
          mscol = set_var_col_const(maskpaths.vars);
      end
  end
end

% check for required vars, quit if they can't be found
check_vars = {'sinfo',{'epidata','rsqrpaths','maskpaths'}};
check_required_vars;

cols = {};

% init outdata
if exist('rsqrpaths','var')
  cols = [cols set_var_col_const(rsqrpaths.vars)];
  % rsqr files
  outdata.vars = [outdata.vars 'rsqrpaths'];
  ridx = length(outdata.vars);
  outdata.data{ridx} = rsqrpaths;
end

if exist('residata','var')
  cols = [cols set_var_col_const(residata.vars)];
  % residuals
  outdata.vars = [outdata.vars 'resid'];
  res_idx = length(outdata.vars);
  outdata.data{res_idx} = residata;
end

if exist('epidata','var')
  cols = [cols set_var_col_const(epidata.vars)];
  % epidata
  outdata.vars = [outdata.vars 'epi'];
  epi_idx = length(outdata.vars);
  outdata.data{epi_idx} = epidata;
end

% subject info
outdata.vars = [outdata.vars 'sinfo'];
sinfo_idx = length(outdata.vars);
outdata.data{sinfo_idx} = ensemble_init_data_struct();
outdata.data{sinfo_idx}.type = 'sinfo';
outdata.data{sinfo_idx}.data = sinfo;

% 
% ? transform images from subject space to group space ?
% 

coreg_list = defs.fmri.general.coreg_list;

for ij=1:length(outdata.data)-1
  group_img_list = outdata.data{ij}.data{cols{ij}.path};
  subids = outdata.data{ij}.data{cols{ij}.subject_id};
  for ii=1:length(group_img_list)
    lfile = group_img_list{ii};
    
    fprintf(1,'processing: %s\n',lfile);
    [fp,fn,fx] = fileparts(lfile);
    nfile = fullfile(fp,sprintf('%s_xfm%s',fn,fx));
  
    lsub  = subids{ii};
  
    xfm_path = fullfile(defs.paths.inroot,lsub,'session1','xfm');
  
    c2t_path = fullfile(xfm_path,sprintf('%s_coplanar2template.mat',lsub));
    if ~exist(c2t_path,'file')
      fprintf(1,'generating cop 2 temp xfm\n');
      h2t_path = fullfile(xfm_path,sprintf('%s_hires2template.mat',lsub));
      c2h_path = fullfile(xfm_path,sprintf('%s_coplanar2hires.mat',lsub));
      xfmstr = sprintf('convert_xfm -omat %s -concat %s %s',...
          c2t_path,h2t_path,c2h_path);
      status = unix(xfmstr);
      if status, error('problem calculating transformation files'); end
    end
    
    e2t_path = fullfile(xfm_path,sprintf('%s_epi2template.mat',lsub));
    if ~exist(e2t_path,'file')
      fprintf(1,'generating epi 2 temp xfm\n');
      e2c_path = fullfile(xfm_path,sprintf('%s_epi2coplanar.mat',lsub));
      xfmstr = sprintf('convert_xfm -omat %s -concat %s %s',...
          e2t_path,c2t_path,e2c_path);
      status = unix(xfmstr);
      if status, error('problem calculating transformation files'); end
    end
    
    runstridx = strfind(lfile,'run');
    if ~isempty(runstridx)
      runnum = str2num(lfile(runstridx(end)+3));
      if runnum > 1
        %%%% NOTE: will croak for studies with more than 9 runs
        %%%% NOTE: assumes reference run is the first run
        fprintf(1,'generating run %d 2 temp xfm\n',runnum);
        r2t_path = fullfile(xfm_path,sprintf('%s_run%d_2_template.mat',...
            lsub,runnum));
        if ~exist(r2t_path,'file')
          r2r_path = fullfile(xfm_path,sprintf('%s_run%d_run%d.mat',...
              lsub,runnum,1));
          xfmstr = sprintf('convert_xfm -omat %s -concat %s %s',...
              r2t_path,e2t_path,r2r_path);
          status = unix(xfmstr);
          if status, error('problem calculating transformation files'); end
        end
        
        tgt_path = r2t_path;
      else
        tgt_path = e2t_path;
      end
    else
      tgt_path = e2t_path;
    end
    fprintf(1,'applying %s transformation to %s\n',tgt_path,lfile);
    fstr = sprintf('flirt -in %s -ref %s -applyxfm -init %s -out %s',...
        lfile,defs.fmri.fsl.paths.hires_template_path,tgt_path,nfile);
    status = unix(fstr);
    if status, error('problem calculating transformation files'); end
    outdata.data{ij}.data{cols{ij}.path}{ii} = nfile;
  end
end % for ij=1:length(outdata.data{
