function outdata = ensemble_fmri_init(indata,defs)

% inits an fmri dataset on disk for indata.data{inc.subjec_id}, returns paths to relevant files
% 
% REQUIRES
%   defs.sinfo(:)
%   defs.subids
%   defs.paths.inroot
%   defs.paths.outroot
%   defs.expinfo.id
%   defs.fmri.protocol.id
%   defs.init.CHECK_HEADERS (def: 1)
%   defs.init.CORRECT_HIRES (def: 1)
%   defs.init.SYMLINK_COPLANARS (def: 1)
%   defs.init.SYMLINK_EPIS (def: 1)
%   defs.init.VOL_CHECK (def: 1)
%   defs.init.USER_SCANNER_MOCO
%   defs.init.MAKE_4D_NIFTI
%   defs.init.USE_4D_NIFTI
%   defs.init.TOUCH_HEADERS
%   defs.init.REPLACE_BAD_VOLS
%   defs.init.ROTATE_EPI
%   defs.init.WRITE2FILE
%       if this, and defs.init.VOL_CHECK are set, then you must include
%       defs.paths.figpath
%   defs.init.USE_SPM
%   defs.init.CLOBBER
% 
% 2008/08/11 FB - adapted from proc_nostalgia_fmri_fmri

outdata = ensemble_init_data_struct();

global r
VERBOSE = 1;

r = init_results_struct;

r.type = 'fmri_anal';  % Identify the type of this reporting instance
r.report_on_fly = 1;

%%% INITIALIZE REPORTING VAR/FUNC

% get subids
if isfield(defs,'subids')
  proc_subs = defs.subids;
  if ~ischar(proc_subs) && ~iscell(proc_subs)
    fprintf(1,'please provide valid subids');
    return
  elseif ~iscell(proc_subs)
    proc_subs = {proc_subs};
  end
else
  fprintf(1,'please provide valid subids');
  return
end

exp_inroot = defs.paths.inroot;

% Check to make sure that the destination directory exists and create it if necessary
exp_outroot = defs.paths.outroot;
check_dir(exp_outroot);

try CHECK_HEADERS = defs.init.CHECK_HEADERS;
catch CHECK_HEADERS = 1; end
try CORRECT_HIRES = defs.init.CORRECT_HIRES;
catch CORRECT_HIRES = 1; end
try SYMLINK_COPLANARS = defs.init.SYMLINK_COPLANARS;
catch SYMLINK_COPLANARS = 1; end
try SYMLINK_EPIS = defs.init.SYMLINK_EPIS;
catch SYMLINK_EPIS = 1; end
try VOL_CHECK = defs.init.VOL_CHECK;
catch VOL_CHECK = 1; end
try USE_SCANNER_MOCO = defs.init.USE_SCANNER_MOCO;
catch USE_SCANNER_MOCO = 0; end
try MAKE_4D_NIFTI = defs.init.MAKE_4D_NIFTI;
catch MAKE_4D_NIFTI = 0; end
try USE_4D_NIFTI = defs.init.USE_4D_NIFTI;
catch USE_4D_NIFTI = 0; end
try TOUCH_HEADERS = defs.init.TOUCH_HEADERS;
catch TOUCH_HEADERS = 0; end
try REPLACE_BAD_VOLS = defs.init.REPLACE_BAD_VOLS;
catch REPLACE_BAD_VOLS = 0; end
try ROTATE_EPI = defs.init.ROTATE_EPI;
catch ROTATE_EPI = 0; end
try USE_SPM = defs.init.USE_SPM; catch USE_SPM = 0; end
try CLOBBER = defs.init.CLOBBER; catch CLOBBER = 0; end
try WRITE2FILE = defs.init.WRITE2FILE; catch WRITE2FILE = 0; end

% set up data structs
if isfield(defs,'sinfo') && isstruct(defs.sinfo)
  sinfo = defs.sinfo;
end

% check for required vars
check_vars = {'sinfo'};
check_required_vars;

outdata.vars = [outdata.vars 'sinfo'];
sinfo_idx = length(outdata.vars);
outdata.data{sinfo_idx} = ensemble_init_data_struct();
outdata.data{sinfo_idx}.type = 'sinfo';
outdata.data{sinfo_idx}.data = sinfo;

outdata.vars = [outdata.vars 'paths'];
paths_idx = length(outdata.vars);
outdata.data{paths_idx} = ensemble_init_data_struct();
outdata.data{paths_idx}.type = 'paths';
outdata.data{paths_idx}.vars = {'subject_id','session','run',...
    'path_type','path'};
outdata.data{paths_idx}.data{1} = {};
outdata.data{paths_idx}.data{2} = [];
outdata.data{paths_idx}.data{3} = [];
outdata.data{paths_idx}.data{4} = {};
outdata.data{paths_idx}.data{5} = {};
pathcol = set_var_col_const(outdata.data{paths_idx}.vars);

if CORRECT_HIRES
  outdata.vars = [outdata.vars 'hires'];
  hires_idx = length(outdata.vars);
  outdata.data{hires_idx} = ensemble_init_data_struct();
  outdata.data{hires_idx}.type='hires';
  outdata.data{hires_idx}.vars = {'subject_id','session',...
      'ensemble_id','path'};
  outdata.data{hires_idx}.data{1} = {};
  outdata.data{hires_idx}.data{2} = [];
  outdata.data{hires_idx}.data{3} = [];
  outdata.data{hires_idx}.data{4} = {};
  hscol = set_var_col_const(outdata.data{hires_idx}.vars);
end

if SYMLINK_COPLANARS
  outdata.vars = [outdata.vars 'coplanar'];
  cop_idx = length(outdata.vars);
  outdata.data{cop_idx} = ensemble_init_data_struct();
  outdata.data{cop_idx}.type='coplanar';
  outdata.data{cop_idx}.vars = {'subject_id','session',...
      'ensemble_id','path'};
  outdata.data{cop_idx}.data{1} = {};
  outdata.data{cop_idx}.data{2} = [];
  outdata.data{cop_idx}.data{3} = [];
  outdata.data{cop_idx}.data{4} = {};
  cocol = set_var_col_const(outdata.data{cop_idx}.vars);
end

if SYMLINK_EPIS || MAKE_4D_NIFTI
  outdata.vars = [outdata.vars 'epi'];
  epi_idx = length(outdata.vars);
  outdata.data{epi_idx} = ensemble_init_data_struct();
  outdata.data{epi_idx}.type='epi';
  outdata.data{epi_idx}.vars = {'subject_id','session',...
      'ensemble_id','run','path'};
  outdata.data{epi_idx}.data{1} = {};
  outdata.data{epi_idx}.data{2} = [];
  outdata.data{epi_idx}.data{3} = [];
  outdata.data{epi_idx}.data{4} = [];
  outdata.data{epi_idx}.data{5} = {};
  epcol = set_var_col_const(outdata.data{epi_idx}.vars);
end

outcols = set_var_col_const(outdata.vars);

%
% START OF THE SUBJECT LOOP
%

nsub_proc = length(proc_subs);

for isub=1:nsub_proc
  subidx = strmatch(proc_subs{isub},{sinfo.id},'exact');
  subid = sinfo(subidx).id;
  msg = sprintf('\t\tPROCESSING SUBJECT (%d/%d): %s\n', isub, nsub_proc,subid);
  r = update_report(r,msg);
  
  % Deal with directory infrastructure
  sub_indir = fullfile(exp_inroot, subid);
  sub_outdir = fullfile(exp_outroot, subid);
  check_dir(sub_outdir,1);

  outdata.data{paths_idx}.data{1} = [outdata.data{paths_idx}.data{1};...
      subid];
  outdata.data{paths_idx}.data{2} = [outdata.data{paths_idx}.data{2};...
      0];
  outdata.data{paths_idx}.data{3} = [outdata.data{paths_idx}.data{3};...
      0];
  outdata.data{paths_idx}.data{4} = [outdata.data{paths_idx}.data{4};...
      'sub_indir'];
  outdata.data{paths_idx}.data{5} = [outdata.data{paths_idx}.data{5};...
      sub_indir];

  outdata.data{paths_idx}.data{1} = [outdata.data{paths_idx}.data{1};...
      subid];
  outdata.data{paths_idx}.data{2} = [outdata.data{paths_idx}.data{2};...
      0];
  outdata.data{paths_idx}.data{3} = [outdata.data{paths_idx}.data{3};...
      0];
  outdata.data{paths_idx}.data{4} = [outdata.data{paths_idx}.data{4};...
      'sub_outdir'];
  outdata.data{paths_idx}.data{5} = [outdata.data{paths_idx}.data{5};...
      sub_outdir];

  % Determine number of sessions for this subject
  nsess = length(sinfo(subidx).sessinfo);
  
  %
  % START OF THE SESSION LOOP
  %
  
  for isess = 1:nsess
    sess = sinfo(subidx).sessinfo(isess);
    
    if ~sess.use_session
      msg = sprintf('\t\t\tSkipping session %d\n', isess);
      r = update_report(r,msg);
      continue
    end
    
    session_stub = sess.id;
    r = update_report(r,sprintf('\t\t\t%s\n', session_stub));
    
    % Determine how many runs we're dealing with
    nruns = length(sess.use_epi_runs);  
    ngood_runs = nruns;
    
    % Deal with the path definitions
    fmri_sess_paths;
    
    % Figure out which version of the experiment applies
    exp_id = sess.exp_id;
    expidx = strmatch(exp_id,{defs.expinfo.id},'exact');
    if isempty(expidx)
      msg = sprintf('!!! Could not match experiment ID: %s\n', sess.exp_id);
      r = update_report(r,msg);
      continue
    end
    
    % Get the scanner protocol for this session
    protocol_idx = ...
	   strmatch(sess.protocol_id,{defs.fmri.protocol.id},'exact');
    protocol = defs.fmri.protocol(protocol_idx);
            
    %
    % Deal with the physiological files for this session
    %
    
    % Here we need to set things up that we need to accomplish at the
    % session level, i.e. across runs
    
    % Swap dimensions on hires images if necessary and if a swap sequence
    % is specified
    if CORRECT_HIRES & ~isempty(protocol.hires.swapseq)
      fmri_correct_hires;
    end % if CORRECT_HIRES
    
    if SYMLINK_COPLANARS
      fmri_symlink_coplanars;
    end % if SYMLINK_COPLANARS
	    
    % Figure out which EPI series we are dealing with
    if USE_SCANNER_MOCO
      series_map_idx = find(ismember([sess.series_mappings{:,2}],{'epi_moco','epi_dico_moco'}));
    else
      series_map_idx = find(ismember([sess.series_mappings{:,2}],{'epi','epi_dico'}));
    end

    if isempty(series_map_idx)
      msg = sprintf('Did not find EPI directory mapping\n');
      r = update_report(r,msg);
      continue
    elseif length(series_map_idx) > 1
      msg = sprintf(['Found multiple EPI series mappings (%s). Please be more' ...
	    ' specific\n'], ...
	  cell2str([sess.series_mappings{series_map_idx,2}],','));
      r = update_report(r,msg);
      continue;
    end

    %
    % START OF RUN LOOP
    %
    nvols = [];
    for irun = 1:nruns
      % Deal with directory naming
      % Run stubs will vary because in the original directories, they have
      % series numbers attached.
      runidx = sess.use_epi_runs(irun);
      
      if ~USE_SCANNER_MOCO
        run_stub = sprintf('%s_%s',...
            sess.series_mappings{series_map_idx,1}{runidx},...
            protocol.epi.dirstub);
      else
        run_stub = '';
        msg = sprintf('Have not handled this option yet\n');
        r = update_report(r,msg);
        continue
      end
      
      run_indir = fullfile(epi_indir, run_stub); % location of original data
      
      if (isempty(dir(run_indir)))
        run_indir = [run_indir '*'];
        run_stub = dir(run_indir);
        if isempty(run_stub)
            msg = sprintf('couldn''t find run_stub (%s)',run_indir);
            r = update_report(r,msg);
            continue
        else
            run_indir = fullfile(epi_indir,run_stub.name); % new run_stub
        end
      end

      run_outdir = fullfile(epi_outdir, sprintf('%s_run%d',subid,irun)); % location of modified data

      % Make sure the output directory exists
      check_dir(run_outdir,1);
        
      outdata.data{paths_idx}.data{1} = [outdata.data{paths_idx}.data{1};...
        subid];
      outdata.data{paths_idx}.data{2} = [outdata.data{paths_idx}.data{2};...
        isess];
      outdata.data{paths_idx}.data{3} = [outdata.data{paths_idx}.data{3};...
        irun];
      outdata.data{paths_idx}.data{4} = [outdata.data{paths_idx}.data{4};...
        'run_indir'];
      outdata.data{paths_idx}.data{5} = [outdata.data{paths_idx}.data{5};...
        run_indir];

      outdata.data{paths_idx}.data{1} = [outdata.data{paths_idx}.data{1};...
        subid];
      outdata.data{paths_idx}.data{2} = [outdata.data{paths_idx}.data{2};...
        isess];
      outdata.data{paths_idx}.data{3} = [outdata.data{paths_idx}.data{3};...
        irun];
      outdata.data{paths_idx}.data{4} = [outdata.data{paths_idx}.data{4};...
        'run_outdir'];
      outdata.data{paths_idx}.data{5} = [outdata.data{paths_idx}.data{5};...
        run_outdir];

      epifstubfmt = protocol.epi.fstubfmt;
            
      % Only execute analyses if this run exists or if we are dealing with a
      % model that is using residuals from a previous model, rather than the
      % EPI data
      try using_resid = strcmp(curr_model.srcdata,'residuals'); catch ...
	    using_resid = 0; end
      
      if ~isempty(dir(fullfile(run_outdir,'*.img'))) || ...
	    ~isempty(dir(fullfile(run_indir,'*.img'))) || using_resid
	

	
	    if SYMLINK_EPIS
          fmri_symlink_epis;
	    end % if SYMLINK_EPIS

	    % Check for consistency in ANALYZE file headers
	    if CHECK_HEADERS
	      fmri_check_headers;
	    end % if CHECK_HEADERS

	    if MAKE_4D_NIFTI
	      fmri_make_4D_nifti;
	    end % if MAKE_4D_NIFTI

  	    if VOL_CHECK
	      fmri_vol_check;
	    end % if VOL_CHECK

	    % Get file lists
	    if USE_4D_NIFTI
	      msg = sprintf('Using 4D Nifti file\n');
	      r = update_report(r,msg);
	      flist = fullfile(run_outdir,sprintf('%s_run%d.nii', subid, irun));
        else
	      %msg = sprintf('Using ANALYZE format files\n');
	      srcdir = run_outdir;  % need to use this if doing anything with
	      % SPM
	      srcstub = sprintf('%s*.img',subid);
	      flist = get_spm_flist(srcdir,srcstub);
  	    end % if USE_4D_NIFTI
	
	    if TOUCH_HEADERS
	      fmri_touch_headers;
	    end % if TOUCH_HEADERS
	
	    if REPLACE_BAD_VOLS & ~isempty(sinfo(subidx).badvols{irun})
	      fmri_replace_bad_vols;
	    end % if REPLACE_BAD_VOLS
	
	    if ROTATE_EPI
	      fmri_rotate_epi;
	    end % if ROTATE_EPI

        if ~exist('flist') || isempty(flist)
          if ~isempty(dir(run_outdir))
            flist = dir(run_outdir);
          elseif ~isempty(dir(run_indir))
            flist = dir(run_indir);
          end
        end
      
        if exist('flist') && ~isempty(flist)
          for ifl=1:length(flist)
            outdata.data{epi_idx}.data{epcol.subject_id} = ...
              [outdata.data{epi_idx}.data{epcol.subject_id}; subid];
            outdata.data{epi_idx}.data{epcol.session} = ...
              [outdata.data{epi_idx}.data{epcol.session}; isess];
            outdata.data{epi_idx}.data{epcol.ensemble_id} = ...
              [outdata.data{epi_idx}.data{epcol.ensemble_id}; sess.ensemble_id];
            outdata.data{epi_idx}.data{epcol.run} = ...
              [outdata.data{epi_idx}.data{epcol.run}; irun];
            outdata.data{epi_idx}.data{epcol.path} = ...
              [outdata.data{epi_idx}.data{epcol.path}; flist(ifl,:)];
          end
        end

      end % if exist_epi
    end % for irun
  end % for isess
end % for isub=
