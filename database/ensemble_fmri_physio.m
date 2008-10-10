function outdata = ensemble_fmri_physio(indata,defs)

% processes physiological data, stores .mat file on disk in run/regressors
% 
% REQUIRES
%   defs.sinfo(:).sessinfo(:).physio
%   defs.sinfo(:).sessinfo(:).physio.source
%   defs.init.WRITE2FILE
% 
% 2008/10/06 FB - started coding

outdata = ensemble_init_data_struct();

global r

r = init_results_struct;

r.type = 'fmri_physio';  % Identify the type of this reporting instance
r.report_on_fly = 1;

%%% INITIALIZE VARIABLES
try WRITE2FILE = defs.init.WRITE2FILE; catch WRITE2FILE = 1; end

% Parse out the input data
for idata = 1:length(indata)
  switch indata(idata).type
    case {'sinfo'}
      sinfo = indata(idata);
      proc_subs = defs.sinfo(:).id;
    case {'paths'}
      pathdata = indata(idata);
      pcol = set_var_col_const(pathdata.vars);
  end
end

% check for required vars
check_vars = {'sinfo','pathdata'};
check_required_vars;

outdata.vars = [outdata.vars 'sinfo'];
sinfo_idx = length(outdata.vars);
outdata.data{sinfo_idx} = ensemble_init_data_struct();
outdata.data{sinfo_idx}.type = 'sinfo';
outdata.data{sinfo_idx}.data = sinfo;

outdata.vars = [outdata.vars 'physio_data'];
pdata_idx = length(outdata.vars);
outdata.data{pdata_idx} = ensemble_init_data_struct();
outdata.data{pdata_idx}.type = 'physio_data';
outdata.data{pdata_idx}.vars = {'subject_id','session',...
    'ensemble_id','physiodata'};
outdata.data{pdata_idx}.data{1} = {};
outdata.data{pdata_idx}.data{2} = [];
outdata.data{pdata_idx}.data{3} = [];
outdata.data{pdata_idx}.data{4} = ensemble_init_data_struct();
pdcol = set_var_col_const(outdata.data{pdata_idx}.vars);

if WRITE2FILE
  outdata.vars = [outdata.vars 'physio_paths'];
  ppaths_idx = length(outdata.vars);
  outdata.data{ppaths_idx} = ensemble_init_data_struct();
  outdata.data{ppaths_idx}.type='physio_paths';
  outdata.data{ppaths_idx}.vars = {'subject_id','session',...
      'ensemble_id','path'};
  outdata.data{ppaths_idx}.data{1} = {};
  outdata.data{ppaths_idx}.data{2} = [];
  outdata.data{ppaths_idx}.data{3} = [];
  outdata.data{ppaths_idx}.data{4} = {};
  ppcol = set_var_col_const(outdata.data{ppaths_idx}.vars);
end

ecdparams.outDataName = 'physio_data';

%
% START OF THE SUBJECT LOOP
%

nsub_proc = length(proc_subs);

for isub=1:nsub_proc
  subidx = strmatch(proc_subs{isub},{sinfo.id},'exact');
  subid = sinfo(subidx).id;
  msg = sprintf('\t\tPROCESSING SUBJECT (%d/%d): %s\n', isub, nsub_proc,subid);
  r = update_report(r,msg);

  % Determine number of sessions for this subject
  nsess = length(sinfo(subidx).sessinfo);
  
  %
  % START OF THE SESSION LOOP
  %
  
  for isess = 1:nsess
    sess = sinfo(subidx).sessinfo(isess);
    
    if isfield(sess,'use_session') && ~sess.use_session
      msg = sprintf('\t\t\tSkipping session %d\n', isess);
      r = update_report(r,msg);
      continue
    elseif ~isfield(sess,'physio') || ~isstruct(sess.physio)
      msg = sprintf(['\t\t\tNo Physio sinfo for session %d, '...
          'SKIPPING!\n'],isess);
      r = update_report(r,msg);
      continue
    end
    
    session_stub = sess.id;
    r = update_report(r,sprintf('\t\t\t%s\n', session_stub));
    
    % init session vars
    sdata = ensemble_init_data_struct();
    physio = sinfo(subidx).sessinfo(isess).physio;
    pparams = init_physmon_params(physio.source);
    pparams.ftypes = {physio.channels(:).name};
    
    % Obtain information about the scanning protocol that was used in this
    % session
    protocol_idx = strmatch(sess.protocol_id, {defs.fmri.protocol.id}, ...
	'exact');
    
    pparams.TR = defs.fmri.protocol(protocol_idx).epi.tr;
    pparams.nslice_per_vol = defs.fmri.protocol(protocol_idx).epi.nslices;
    
    % get physio_indir/outdir
    sfilt = struct();
    sfilt.include.all.subject_id = {subid};
    sfilt.include.all.session = isess;
    spaths = ensemble_filter(pathdata,spaths);
    
    indx = strmatch('physio_indir',spaths.data{pcol.path});
    physio_indir = spaths.data{pcol.path}{indx};

    outdx = strmatch('physio_outdir',spaths.data{pcol.path});
    physio_outdir = spaths.data{pcol.path}{outdx};
    
    % Determine how many runs we're dealing with
    runs = sess.use_epi_runs;
    nruns = length(runs);  
    
    %
    % START OF RUN LOOP
    %
    for irun = 1:nruns

      physiofname = fullfile(physio_indir,physio.datafiles{runs(irun)});

      % % % % % DEAL WITH PHYSIO DATA FILES
      
      % % % % % CONCATENATE PHYSIO DATA WITH SDATA
      
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
      if isempty(outdata.data{pdata_idx}.data{pdcol.physiodata}.vars)
        outdata.data{pdata_idx}.data{pdcol.physiodata}.vars = sdata.vars;
        for iv = 1:length(sdata.vars)
          if ischar(sdata.data{iv})
            outdata.data{pdata_idx}.data{pdcol.physiodata}.data{iv} = {};
          else
            outdata.data{pdata_idx}.data{pdcol.physiodata}.data{iv} = [];
          end
        end
        outdata.data{pdata_idx}.data{pdcol.physiodata}.data = ...
            ensemble_concat_datastruct(...
            {sdata,outdata.data{pdata_idx}.data{pdcol.physiodata}.data},...
            ecdparams);
      end
    end
    
    % write out sdata
    if WRITE2FILE
      mat_fname = fullfile(physio_outdir,sprintf('%s_sess%d_physio.mat',...
          subid,isess));
      save(mat_fname,'-struct','sdata');
    end
    
  end % for isess
end % for isub=
