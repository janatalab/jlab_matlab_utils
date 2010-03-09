function outdata = ensemble_fmri_vol_check(indata,defs)

% loads epi data, plots mean slice activation across time
% 
%   outdata = ensemble_fmri_vol_check(indata,defs)
% 
% REQUIRES
%   indata
%       epidata
%   defs
% 
% RETURNS
%   outdata
%       epidata
%       figpaths
% 
% 2010.02.16 FB - adapted from private/matlab/projects/...
%   autobio/fmri/fmri_vol_check.m

global r

outdata = ensemble_init_data_struct();
outdata.type = 'fmri_vol_check';
outdata.name = outdata.type;

r = init_results_struct;

r.type = 'fmri_vol_check';  % Identify the type of this reporting instance
r.report_on_fly = 1;

% Parse out the input data
for idata = 1:length(indata)
  if isfield(indata{idata},'type')
    switch indata{idata}.type
      case 'epi'
        epidata = indata{idata};
        epicol = set_var_col_const(epidata.vars);
      case 'paths'
        pathdata = indata{idata};
        pacol = set_var_col_const(pathdata.vars);
    end
  end
end

% check for required vars
check_vars = {'epidata','pathdata'};
check_required_vars;

if (iscell(indata) && ~isempty(indata) && isfield(indata{1},'task') && ...
        ~isempty(strmatch('return_outdir',indata{1}.task))) || ...
        (isstruct(indata) && isfield(indata,'task') && ...
        ~isempty(strmatch('return_outdir',indata.task)))
  if exist('pathdata','var') && ~isempty(pathdata.data{1})
    if length(nsub_proc) == 1
      pfilt = struct();
      pfilt.include.all.subject_id = proc_subs;
      lpathdata = ensemble_filter(pathdata,pfilt);
      if ~isempty(lpathdata.data{1})
        sfilt = pfilt;
        sfilt.include.all.path_type = {'epi_outdir'};
        spathdata = ensemble_filter(lpathdata,sfilt);
        if length(spathdata.data{1}) == 1
          % one epi outdir, save outdata = epi_outdir
          outdata = spathdata.data{pacol.path}{1};
        else
          sfilt = pfilt;
          sfilt.include.all.path_type = {'sess_outdir'};
          spathdata = ensemble_filter(lpathdata,sfilt);
          if length(spathdata.data{1}) == 1;
            outdata = spathdata.data{pacol.path}{1};
          else
            sfilt = pfilt;
            sfilt.include.all.path_type = {'sub_outdir'};
            spathdata = ensemble_filter(lpathdata,sfilt);
            if length(spathdata.data{1}) == 1;
              outdata = spathdata.data{pacol.path}{1};            
            end
          end
        end
      end
    end
  end
  if ~exist('outdata','var') || ~exist(outdata,'dir'), outdata = ''; end
  return
end

% epidata
outdata.vars = [outdata.vars 'epi'];
epi_idx = length(outdata.vars);
outdata.data{epi_idx} = epidata;

% figure path output data structure
outdata.vars = [outdata.vars 'figpaths'];
paths_idx = length(outdata.vars);
outdata.data{paths_idx} = ensemble_init_data_struct();
outdata.data{paths_idx}.type = 'figpaths';
outdata.data{paths_idx}.vars = {'subject_id','session','path'};
fpcol = set_var_col_const(outdata.data{paths_idx}.vars);
outdata.data{paths_idx}.data{fpcol.subject_id} = {};
outdata.data{paths_idx}.data{fpcol.session} = [];
outdata.data{paths_idx}.data{fpcol.path} = {};

% init flags
try WRITE2FILE = defs.vol_check.WRITE2FILE; catch WRITE2FILE = 1; end
try figstub = ['_' defs.vol_check.figstub]; catch figstub = ''; end
try popts = defs.vol_check.printargs; catch popts = {'-dpsc'}; end

% iterate over subjects
subids = unique(epidata.data{epicol.subject_id});
nsub = length(subids);
for j=1:nsub
  subid = subids{j};
  subfilt.include.all.subject_id = {subid};
  subdata = ensemble_filter(epidata,subfilt);
    
  usess = unique(subdata.data{epicol.session});
  nsess = length(usess);

  % iterate over sessions
  for k=1:nsess
    sid = usess(k);
    sfilt.include.all.session = sid;
    sdata = ensemble_filter(subdata,sfilt);
    
    figfname = fullfile(defs.paths.figpath,...
        sprintf('fmri_vol_check_%s_sess%d%s.ps',subid,sid,figstub));
    outdata.data{paths_idx} = ensemble_add_data_struct_row(...
        outdata.data{paths_idx},'subject_id',subid,'session',sid,...
        'path',figfname);

    urun = unique(sdata.data{epicol.run});
    nrun = length(urun);

    % iterate over runs
    for l=1:nrun
      rid = urun(l);
      rfilt.include.all.run = rid;
      rdata = ensemble_filter(sdata,rfilt);
      
      flist = rdata.data{epicol.path};
      npts = length(flist);

      % Get dimensions of data
      warning off
      V = spm_vol(flist{1});
      warning on

      % Presize Y which will hold the EPI data
      Y = zeros(V.dim(1),V.dim(2),V.dim(3),npts);
      for ipt = 1:npts
        warning off
        V = spm_vol(flist{ipt});
        warning on
        msg = sprintf('.');
        r = update_report(r,msg);
        Y(:,:,:,ipt) = spm_read_vols(V);
      end
      msg = sprintf('\n');
      r = update_report(r,msg);      

      msg = sprintf('Calculating mean intensity for run %d\n',rid);
      r = update_report(r,msg);

      % Calculate summed intensity within each slice
      % Method 1 for summing across 2 dimensions
      Ymean = mean(Y,1);
      Ymean = mean(Ymean,2);
      Ymean = squeeze(Ymean);

      subplot(nrun,1,l)
      imagesc(Ymean), colorbar

      set(gca,'xtick',0:fix(npts/10):npts)

      ylabel('Slice#')
      xlabel('Volume#')

      if l == 1
        title(sprintf('Subject %s: Session %d: Run %d',subid,sid,rid));
        lpopts = popts;
      else
        title(sprintf('Run %d',rid));
        lpopts = {popts{:},'-append'};
      end

      if WRITE2FILE
        msg = sprintf('Writing mean intensity plot to %s',figfname);
        r = update_report(r,msg);
        print(figfname,lpopts{:});
      end

    end % for l=1:nrun
  end % for k=1:nsess
end % for j=1:nsub


