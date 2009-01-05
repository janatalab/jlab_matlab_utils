function outdata = ensemble_fmri_parse_present(indata,defs)

% parses presentation files, saves to appropriate place on disk
% 
% REQUIRES
%   defs.sinfo(:).sessinfo(:).pres
%   defs.paths.inroot
%   defs.paths.outroot
%   defs.expinfo.id
%   defs.fmri.protocol.id
%   defs.init.WRITE2FILE
%       if this, and defs.init.VOL_CHECK are set, then you must include
%       defs.paths.figpath
%   defs.init.USE_SPM
%   defs.init.CLOBBER
% 
% 2008/08/11 FB - adapted from proc_nostalgia_fmri_fmri

outdata = ensemble_init_data_struct();

global r

r = init_results_struct;

r.type = 'fmri_anal';  % Identify the type of this reporting instance
r.report_on_fly = 1;

%%% INITIALIZE VARIABLES
try WRITE2FILE = defs.present.WRITE2FILE; catch WRITE2FILE = 1; end

% Parse out the input data
for idata = 1:length(indata)
  if isfield(indata{idata},'type')
    switch indata{idata}.type
      case {'paths'}
        pathdata = indata{idata};
        pcol = set_var_col_const(pathdata.vars);
    end
  end
end

if isfield(defs,'sinfo')
  sinfo = defs.sinfo;
  proc_subs = {sinfo(:).id};
  nsub_proc = length(proc_subs);
end

% check for required vars
check_vars = {'sinfo','pathdata'};
check_required_vars;

if (iscell(indata) && ~isempty(indata) && isfield(indata{1},'task') && ...
        ~isempty(strmatch('return_outdir',indata{1}.task))) || ...
        (isstruct(indata) && isfield(indata,'task') && ...
        ~isempty(strmatch('return_outdir',indata.task)))
  outdata = defs.paths.outroot;
  if length(nsub_proc) == 1
    outdata = fullfile(outdata,proc_subs{1});
    if length(sinfo(1).sessinfo) == 1
      sdir = fullfile(outdata,'session1');
      if exist(sdir,'dir')
        outdata = sdir;
      end
    else
      % multiple sessions, save in the current outdir
    end
  else
    % multiple subjects, save in defs.paths.outroot
  end
  if ~exist(outdata,'dir'), outdata = ''; end
  return
end

outdata.vars = [outdata.vars 'sinfo'];
sinfo_idx = length(outdata.vars);
outdata.data{sinfo_idx} = ensemble_init_data_struct();
outdata.data{sinfo_idx}.type = 'sinfo';
outdata.data{sinfo_idx}.data = sinfo;

outdata.vars = [outdata.vars 'pres_data'];
pdata_idx = length(outdata.vars);
outdata.data{pdata_idx} = ensemble_init_data_struct();
outdata.data{pdata_idx}.type = 'presentation_data';
outdata.data{pdata_idx}.vars = {'subject_id','session',...
    'ensemble_id','presdata'};
outdata.data{pdata_idx}.data{1} = {};
outdata.data{pdata_idx}.data{2} = [];
outdata.data{pdata_idx}.data{3} = [];
outdata.data{pdata_idx}.data{4} = ensemble_init_data_struct();
pdcol = set_var_col_const(outdata.data{pdata_idx}.vars);

if WRITE2FILE
  outdata.vars = [outdata.vars 'pres_paths'];
  ppaths_idx = length(outdata.vars);
  outdata.data{ppaths_idx} = ensemble_init_data_struct();
  outdata.data{ppaths_idx}.type='pres_paths';
  outdata.data{ppaths_idx}.vars = {'subject_id','session',...
      'ensemble_id','path'};
  outdata.data{ppaths_idx}.data{1} = {};
  outdata.data{ppaths_idx}.data{2} = [];
  outdata.data{ppaths_idx}.data{3} = [];
  outdata.data{ppaths_idx}.data{4} = {};
  ppcol = set_var_col_const(outdata.data{ppaths_idx}.vars);
end

ecdparams.outDataName = 'presentation_data';

%
% START OF THE SUBJECT LOOP
%

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
    
    if isfield(sess,'use_session') && ~sess.use_session
      msg = sprintf('\t\t\tSkipping session %d\n', isess);
      r = update_report(r,msg);
      continue
    elseif ~isfield(sess,'pres') || ~isstruct(sess.pres)
      msg = sprintf(['\t\t\tNo Presentation sinfo for session %d, '...
          'SKIPPING!\n'],isess);
      r = update_report(r,msg);
      continue
    end
    
    session_stub = sess.id;
    r = update_report(r,sprintf('\t\t\t%s\n', session_stub));
    
    % init session vars
    sdata = ensemble_init_data_struct();    
    pres = sinfo(isub).sessinfo(isess).pres;
    
    % get behav_indir/outdir
    sfilt = struct();
    sfilt.include.all.subject_id = {subid};
    sfilt.include.all.session = isess;
    spaths = ensemble_filter(pathdata,sfilt);
    
    indx = strmatch('behav_indir',spaths.data{pcol.path_type});
    behav_indir = spaths.data{pcol.path}{indx};

    outdx = strmatch('behav_outdir',spaths.data{pcol.path_type});
    behav_outdir = spaths.data{pcol.path}{outdx};
    
    andx = strmatch('anal_outdir',spaths.data{pcol.path_type});
    anal_outdir = spaths.data{pcol.path}{andx};
    
    % Determine how many runs we're dealing with
    runs = sess.use_epi_runs;
    nruns = length(runs);  
    
    %
    % START OF RUN LOOP
    %
    for irun = 1:nruns

      presfname = fullfile(behav_indir,pres.logfiles{runs(irun),1});
      targruns  = pres.logfiles{runs(irun),2};

      pdata = presentation_parse(presfname,pres.params);
      if isempty(pdata)
        msg = sprintf('No Presentation Data Returned, SKIPPING\n');
        r = update_report(r,msg);
        continue
      end
      
      PL = set_var_col_const(pdata.vars);
      
      % filter for target runs within the pres file, identified in sinfo
      pfilt = struct();
      pfilt.include.all.RUN = targruns;
      pdata = ensemble_filter(pdata,pfilt);
      
      % set pdata.data{PL.RUN} to the run # in sinfo
      pdata.data{PL.RUN} = ones(length(pdata.data{1}),1)*runs(irun);
      
      % init sdata vars with pdata vars if no info in sdata
      if isempty(sdata.vars)
        sdata.vars = pdata.vars;
        for iv = 1:length(sdata.vars)
          if ischar(pdata.data{iv})
            sdata.data{iv} = {};
          else
            sdata.data{iv} = [];
          end
        end
      end
      
      sdata = ensemble_concat_datastruct({pdata,sdata},ecdparams);
      
    end % for irun
    
    % concatenate sdata with outdata
    if ~isempty(sdata.data{1})
      outdata.data{pdata_idx}.data{pdcol.subject_id} = [...
          outdata.data{pdata_idx}.data{pdcol.subject_id}; subid];
      outdata.data{pdata_idx}.data{pdcol.session} = [...
          outdata.data{pdata_idx}.data{pdcol.session}; isess];
      outdata.data{pdata_idx}.data{pdcol.ensemble_id} = [...
          outdata.data{pdata_idx}.data{pdcol.ensemble_id}; sess.ensemble_id];
      if isempty(outdata.data{pdata_idx}.data{pdcol.presdata}.vars)
        outdata.data{pdata_idx}.data{pdcol.presdata}.vars = sdata.vars;
        for iv = 1:length(sdata.vars)
          if ischar(sdata.data{iv})
            outdata.data{pdata_idx}.data{pdcol.presdata}.data{iv} = {};
          else
            outdata.data{pdata_idx}.data{pdcol.presdata}.data{iv} = [];
          end
        end
        outdata.data{pdata_idx}.data{pdcol.presdata}.data = ...
            ensemble_concat_datastruct(...
            {sdata,outdata.data{pdata_idx}.data{pdcol.presdata}},...
            ecdparams);
      end
    end
    
    % write out sdata
    if WRITE2FILE
      mat_fname = fullfile(behav_outdir,sprintf('%s_sess%d_present.mat',...
          subid,isess));
      save(mat_fname,'-struct','sdata');
      outdata.data{ppaths_idx}.data{ppcol.subject_id} = [...
          outdata.data{ppaths_idx}.data{ppcol.subject_id}; subid];
      outdata.data{ppaths_idx}.data{ppcol.session} = [...
          outdata.data{ppaths_idx}.data{ppcol.session}; isess];
      outdata.data{ppaths_idx}.data{ppcol.ensemble_id} = [...
          outdata.data{ppaths_idx}.data{ppcol.ensemble_id}; sess.ensemble_id];
      outdata.data{ppaths_idx}.data{ppcol.path} = [...
          outdata.data{ppaths_idx}.data{ppcol.path}; mat_fname];
      
      xdefs = defs.present.export_params;
      xdefs.export.fname = fullfile(anal_outdir,...
          sprintf('%s_sess%d_present.csv',subid,isess));
      xdefs.sas.fname = fullfile(anal_outdir,...
          sprintf('%s_sess%d_present.sas',subid,isess));
      
      ensemble_export_sastxt(sdata,xdefs);
    end
    
  end % for isess
end % for isub=
