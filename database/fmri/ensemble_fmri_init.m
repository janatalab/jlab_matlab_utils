function outdata = ensemble_fmri_init(indata,defs)

% inits an fmri dataset on disk for indata.data{inc.subjec_id}, returns paths to relevant files
% 
% REQUIRES
%   defs.sinfo(:) - containing subject information for all subjects you
%       want to process
%   defs.paths.inroot
%   defs.paths.outroot
%   defs.expinfo.id
%   defs.fmri.protocol.id
%   defs.init.CHECK_HEADERS (default: 1)
%   defs.init.CORRECT_HIRES (default: 1)
%   defs.init.SYMLINK_COPLANARS (default: 1)
%   defs.init.SYMLINK_EPIS (default: 1)
%   defs.init.VOL_CHECK (default: 1)
%   defs.init.USER_SCANNER_MOCO
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

% get sinfo
if isfield(defs,'sinfo') && isstruct(defs.sinfo)
  sinfo = defs.sinfo;
  proc_subs = {sinfo(:).id};
  nsub_proc = length(proc_subs);
else
  warning(1,'please provide a valid sinfo struct');
  return
end

exp_inroot = defs.paths.inroot;

% Check to make sure that the destination directory exists and create it if necessary
exp_outroot = defs.paths.outroot;
check_dir(exp_outroot);

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
try TOUCH_HEADERS = defs.init.TOUCH_HEADERS;
catch TOUCH_HEADERS = 0; end
try REPLACE_BAD_VOLS = defs.init.REPLACE_BAD_VOLS;
catch REPLACE_BAD_VOLS = 0; end
try ROTATE_EPI = defs.init.ROTATE_EPI;
catch ROTATE_EPI = 0; end
try USE_SPM = defs.init.USE_SPM; catch USE_SPM = 0; end
try CLOBBER = defs.init.CLOBBER; catch CLOBBER = 1; end
try WRITE2FILE = defs.init.WRITE2FILE; catch WRITE2FILE = 0; end

try CHECK_HEADERS = defs.init.CHECK_HEADERS;
catch
    if SYMLINK_EPIS
        CHECK_HEADERS = 1;
    else
        CHECK_HEADERS = 0;
    end
end


% check for required vars
check_vars = {'sinfo'};
check_required_vars;

% return the default output directory, if ensemble_jobman_parallel asks
if (iscell(indata) && ~isempty(indata) && isfield(indata{1},'task') && ...
        ~isempty(strmatch('return_outdir',indata{1}.task))) || ...
        (isstruct(indata) && isfield(indata,'task') && ...
        ~isempty(strmatch('return_outdir',indata.task)))
  outdata = exp_outroot;
  if length(nsub_proc) == 1
    outdata = fullfile(outdata,proc_subs{1});
  end
  return
end

% 
% initialize output data structures
% 

% % subject info output data structure
% holds an sinfo structure for all subjects processed by this instance of
% this job
outdata.vars = [outdata.vars 'sinfo'];
sinfo_idx = length(outdata.vars);
outdata.data{sinfo_idx} = ensemble_init_data_struct();
outdata.data{sinfo_idx}.type = 'sinfo';
outdata.data{sinfo_idx}.data = sinfo;

% % path output data structure
% holds all paths created by this instance of this job
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
  % % hires output data structure
  % holds paths to all hires images created by this instance of this job
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
  % % coplanar output data structure
  % holds paths to all coplanar images symlinked by this instance of this job
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

if SYMLINK_EPIS
  % % epi output data structure
  % holds paths to all epi images created or symlinked by this instance of
  % this job
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

for isub=1:nsub_proc
  subid = sinfo(isub).id;
  msg = sprintf('\t\tPROCESSING SUBJECT (%d/%d): %s\n',isub,nsub_proc,subid);
  r = update_report(r,msg);
  
  % Deal with directory infrastructure
  sub_indir = fullfile(exp_inroot, subid);
  sub_outdir = fullfile(exp_outroot, subid);
  check_dir(sub_outdir,1);

  % save paths to the path output data structure
  % subject source file directory
  outdata.data{paths_idx} = ensemble_add_data_struct_row(...
      outdata.data{paths_idx},'subject_id',subid,'session',0,'run',0,...
      'path_type','sub_indir','path',sub_indir);

  % subject destination file directory (often the same as the source)
  outdata.data{paths_idx} = ensemble_add_data_struct_row(...
      outdata.data{paths_idx},'subject_id',subid,'session',0,'run',0,...
      'path_type','sub_indir','path',sub_indir);

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
    
    session_stub = sess.id;
    r = update_report(r,sprintf('\t\t\t%s\n', session_stub));
    
    % Determine how many runs we're dealing with
    nruns = length(sess.use_epi_runs);  
    ngood_runs = nruns;
    
    % Deal with the path definitions
    ensemble_fmri_set_sess_paths;
    
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

    % Swap dimensions on hires images if necessary and if a swap sequence
    % is specified
    if CORRECT_HIRES && ~isempty(protocol.hires.swapseq)
        
        %%%%% CORRECT HIRES
        % swap dimensions
        % concatenate slices, if not a 3d image
        
        msg = sprintf('\t\t\tSwapping dimensions on hires image\n');
        r = update_report(r,msg);

        series_map_idx = strmatch('hires',[sess.series_mappings{:,2}],'exact');
        if isempty(series_map_idx)
          msg = sprintf('Did not find hires for this session\n');
          r = update_report(r,msg);
        else
          hires_indir = fullfile(anat_indir,sprintf('%s_%s', ...
              sess.series_mappings{series_map_idx,1}{1}, ...
              protocol.hires.dirstub));
          infname = fullfile(hires_indir,[protocol.hires.fstub '.img']);

          % Handle the situation where the hires file wasn't necessarily
          % created, but instead I have a bunch of slices
          if ~exist(infname)
            msg = sprintf('Could not find %s.img\n', protocol.hires.fstub);
            r = update_report(r,msg);
            msg = sprintf('Checking for presence of individual slices ...\n');
            r = update_report(r,msg);
            dirlist = dir(fullfile(hires_indir,'0*.img'));
            nfiles = length(dirlist);
            if nfiles
              msg = sprintf('Found %d files.  Concatenating into single volume...\n', nfiles);
              r = update_report(r,msg);
              curr_dir = pwd;
              cd(hires_indir)

              % merge the slices into the volume
              fsl_str = 'fslmerge -a volume 0*.img';
              msg = sprintf('%s\n', fsl_str);
              r = update_report(r,msg);
              unix(fsl_str);

              % convert filetype to NIFTI_PAIR
              fsl_str = 'fslchfiletype NIFTI_PAIR volume';
              msg = sprintf('%s\n', fsl_str);
              r = update_report(r,msg);
              unix(fsl_str);

              cd(curr_dir);
            else
              msg = sprintf('No hires found for %s',subid);
              r = update_report(r,msg);
              error(msg);
            end % if nfiles
          end % if ~exist(infname)

          outfname = fullfile(anat_outdir,sprintf('%s_hires.img', subid));
          if exist(outfname)
            msg = sprintf(['existing hires in anat_outdir for subject %s, '...
                'clobbering']);
            r = update_report(r,msg);
            delete(outfname);
          end
          swap_seq = protocol.hires.swapseq;
          status = fsl_swapdims(infname,swap_seq,outfname, VERBOSE);
          if status
            msg = sprintf('WARNING: Failed to swap dimensions on hires image\n');
            r = update_report(r,msg);
          end

          % add hires info to the hires struct
          outdata.data{hires_idx} = ensemble_add_data_struct_row(...
              outdata.data{hires_idx},'subject_id',subid,'session',...
              isess,'ensemble_id',sess.ensemble_id,'path',outfname);

        end % if isempty(series_map_idx)

        %%%%% END CORRECT HIRES
        
    end % if CORRECT_HIRES
    
    if SYMLINK_COPLANARS
        
        %%%%% SYMLINK_COPLANARS
        % symlink coplanar .img files from orig/
        
        series_map_idx = strmatch('coplanar',[sess.series_mappings{:,2}]);
        ncop = length(series_map_idx);

        for icop = 1:ncop
          coplanar_type = sess.series_mappings{series_map_idx(icop),2}{1};
          switch coplanar_type
            case 'coplanar_T1'
              dirstub = protocol.coplanar_T1.dirstub;
            case 'coplanar_T2'
              dirstub = protocol.coplanar_T2.dirstub;
            otherwise
              dirstub = '';
          end
          coplanar_indir = fullfile(anat_indir,sprintf('%s_%s', ...
              sess.series_mappings{series_map_idx(icop),1}{1}, ...
              dirstub));
          infstub = fullfile(coplanar_indir,protocol.hires.fstub);
          outfstub = fullfile(anat_outdir, sprintf('%s_%s', subid, coplanar_type));

          link_file = sprintf('%s.img', outfstub);
          if exist(link_file,'file')
            if CLOBBER
              fprintf('Removing existing link ...\n');
              unix_str = sprintf('rm %s', link_file);
              unix(unix_str);
            else
              warning('link file %s exists, not clobbering',link_file);
              continue
            end
          end

          unix_str = sprintf('ln -s %s.img %s', infstub, link_file);
          status = unix(unix_str);
          if status
            msg = sprintf('WARNING: Problem linking coplanar\n'); 
            r = update_report(r,msg);
            continue
          end

          unix_str = sprintf('cp %s.hdr %s.hdr', infstub, outfstub);
          status = unix(unix_str);
          if status
            msg = sprintf('WARNING: Problem copying coplanar header\n'); 
            r = update_report(r,msg);
            continue
          else
            % add hires info to the hires struct
            outdata.data{cop_idx} = ensemble_add_data_struct_row(...
                outdata.data{cop_idx},'subject_id',subid,'session',...
                isess,'ensemble_id',...
                sess.ensemble_id,'path',sprintf('%s.img',outfstub));
          end

        end % for icop
        %%%%% END SYMLINK_COPLANARS
        
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
        
      outdata.data{paths_idx} = ensemble_add_data_struct_row(...
          outdata.data{paths_idx},'subject_id',subid,'session',isess,...
          'run',irun,'path_type','run_indir','path',run_indir);

      outdata.data{paths_idx} = ensemble_add_data_struct_row(...
          outdata.data{paths_idx},'subject_id',subid,'session',isess,...
          'run',irun,'path_type','run_outdir','path',run_outdir);

      epifstubfmt = protocol.epi.fstubfmt;
            
      % Only execute analyses if this run exists or if we are dealing with a
      % model that is using residuals from a previous model, rather than the
      % EPI data
      try using_resid = strcmp(curr_model.srcdata,'residuals'); catch ...
	    using_resid = 0; end
      
      if ~isempty(dir(fullfile(run_outdir,'*.img'))) || ...
	    ~isempty(dir(fullfile(run_indir,'*.img'))) || using_resid
	

	
	    if SYMLINK_EPIS
          %%%%% SYMLINK_EPIS
          
            ftypes = {'img','hdr'};
            for itype = 1:length(ftypes)
              flist = dir(fullfile(run_indir,sprintf('*.%s',ftypes{itype})));
              for ifile = 1:length(flist)
                targ_file = fullfile(run_indir,flist(ifile).name);
                outfname = sprintf(epifstubfmt,subid,irun,ifile,ftypes{itype});
                link_file = fullfile(run_outdir, outfname);
                msg = '';
                if exist(link_file)
                  delete(link_file);
                  msg = sprintf('deleted existing file %s\n',link_file);
                end

                switch ftypes{itype}
                  case 'img'
                unix_str = sprintf('ln -s %s %s', targ_file, link_file);
                  case 'hdr'
                unix_str = sprintf('cp %s %s', targ_file, link_file);
                end
                msg = [msg sprintf('%s\n', unix_str)];
                r = update_report(r,msg);
                status = unix(unix_str);

                if status
                  error('Failed to generate symlink')
                end
              end
            end
            
          %%%%% END SYMLINK_EPIS
	    end % if SYMLINK_EPIS

	    % Check for consistency in ANALYZE file headers
	    if CHECK_HEADERS
          %%%%% CHECK_HEADERS

            msg = sprintf('Checking headers in %s\n', run_outdir);
            r = update_report(r,msg);
            srcdir = run_outdir;  % need to use this if doing anything with SPM
            flist = get_spm_flist(srcdir);
            V = spm_vol(flist);
            nfiles = length(V);

            % Get the pixel dimensions for each volume
            msg = sprintf('\tpixdims ...');
            r = update_report(r,msg);
            pixdim_mat = zeros(nfiles,4);
            for ifile = 1:nfiles
              pixdim_mat(ifile,:) = diag(V(ifile).mat)';
            end

            % See if any of the pixdims differ
            diff_mat = diff(pixdim_mat);
            if any(diff_mat(:))
              warning('WARNING: Potentially bad pixel dimension information!\n');
            else
              msg = sprintf(' OK\n');
              r = update_report(r,msg);
            end
        
          %%%%% END CHECK_HEADERS
        end % if CHECK_HEADERS

  	    if VOL_CHECK
          %%%%% VOL_CHECK
          
            % Load the data
            msg = sprintf('Loading ANALYZE format files');
            r = update_report(r,msg);
            flist = dir(fullfile(run_outdir,'*.img'));
            npts = length(flist);
            
            % Get dimensions of data
            fname = fullfile(run_outdir, flist(1).name);
            warning off
            V = spm_vol(fname);
            warning on

            % Presize Y which will hold the EPI data
            Y = zeros(V.dim(1),V.dim(2),V.dim(3),npts);
            for ipt = 1:npts
              fname = fullfile(run_outdir, flist(ipt).name);
              warning off
              V = spm_vol(fname);
              warning on
              msg = sprintf('.');
              r = update_report(r,msg);
              Y(:,:,:,ipt) = spm_read_vols(V);
            end
            msg = sprintf('\n');
            r = update_report(r,msg);

            msg = sprintf('Calculating mean intensity for run %d\n', irun);
            r = update_report(r,msg);

            % Calculate summed intensity within each slice
            % Method 1 for summing across 2 dimensions
            Ymean = mean(Y,1);
            Ymean = mean(Ymean,2);
            Ymean = squeeze(Ymean);

            % Method 2 for summing across 2 dimension
            %Yreshape = reshape(Y,[size(Y,1)*size(Y,2) ...
            %                    size(Y,3) size(Y,4)]);
            %Ymean{irun} = mean(Yreshape);
            %Ymean{irun} = squeeze(Ymean{irun});
            if isfield(defs,'figs') && isfield(defs.figs,'display') ...
                    && defs.figs.display
              if irun == 1
                figure(1), clf
              end

              subplot(ngood_runs,1,irun)
              imagesc(Ymean), colorbar
              colormap(jet)

              % only show volumes 100 through 150
              %        set(gca, 'xlim', [100 250]);
              set(gca,'xtick',0:10:npts)

              ylabel('Slice#')
              xlabel('Volume#')

              if irun == 1
                title(sprintf('Subject %s: Session %d: Run %d',...
                    sinfo(isub).id,isess,runidx));
              else
                title(sprintf('Run %d', runidx));
              end

              if exist('WRITE2FILE','var')
                check_dir(defs.paths.figpath,1);
                fvc_fname = fullfile(defs.paths.figpath,sprintf('fmri_vol_check_%s.ps',...
                  subid));
                msg = sprintf('Writing mean intensity plot to %s',fvc_fname);
                r = update_report(r,msg);
                print(fvc_fname,'-dps');
                msg = sprintf('Converting %s to PDF\n', fvc_fname);
                r = update_report(r,msg);
                pdf_fname = fullfile(defs.paths.figpath,sprintf('fmri_vol_check_%s.pdf',...
                    subid));
                unix_str = sprintf('ps2pdf %s %s', fvc_fname, pdf_fname);
                status = unix(unix_str);
                if status
                  msg = sprintf('Failed to create PDF file from : %s\n', fvc_fname);
                  r = update_report(r,msg);
                end
              end
            end

          %%%%% END VOL_CHECK
	    end % if VOL_CHECK

	    % Get file lists
	    srcdir = run_outdir;  % need to use this if doing anything with
	    % SPM
	    srcstub = sprintf('%s*.img',subid);
	    flist = get_spm_flist(srcdir,srcstub);
	
	    if TOUCH_HEADERS
          %%%%% TOUCH_HEADERS

            nfiles = size(flist,1);
            for ifile = 1:nfiles
              curr_fname = flist(ifile,:);
              msg = sprintf('Writing Nifti header for: %s\n', curr_fname);
              r = update_report(r,msg);
              M = spm_get_space(curr_fname);
              spm_get_space(curr_fname,M);
            end

          %%%%% END TOUCH_HEADERS
        end % if TOUCH_HEADERS
	
	    if REPLACE_BAD_VOLS && ~isempty(sinfo(isub).badvols{irun})
          %%%%% REPLACE_BAD_VOLS
            
            % Make sure we have a directory into which we can copy the
            % original volumes
            orig_dir = fullfile(run_outdir,'orig_badvols');
            check_dir(orig_dir,1);

            nbad = size(sinfo(isub).badvols{irun},1);
            for ibad = 1:nbad
              badvol_idx = sinfo(isub).badvols{irun}{ibad,1};
              goodvol_idxs = sinfo(isub).badvols{irun}{ibad,2};

              badfname = fullfile(run_outdir,sprintf(epifstubfmt, isub, ...
                  badvol_idx));

              % Copy the original bad file
              unix_str = sprintf('cp %s %s', badfname, orig_dir);
              unix(unix_str);

              ngood = length(goodvol_idxs);
              goodflist = {};
              for igood = 1:ngood
                goodflist{igood} = fullfile(run_outdir,sprintf(...
                    epifstubfmt,isub,goodvol_idxs(igood)));
              end

              % Read in the header info for the good data
              warning off
              Vbad = spm_vol(badfname);
              Vgood = spm_vol(strvcat(goodflist));
              warning on
              gooddata = spm_read_vols(Vgood);

              msg = sprintf('Replacing %s with average of volumes %s\n',...
                  badfname,sprintf('%d ', goodvol_idxs));
              r = update_report(r,msg);

              % Remove the old bad file (in case it is a symlink)
              unix_str = sprintf('rm %s', badfname);
              msg = sprintf('%s\n', unix_str);
              r = update_report(r,msg);
              unix(unix_str);

              % Write out average of good images
              Vnew = Vbad;
              Vnew = spm_write_vol(Vnew,mean(gooddata,4)); 

            end % for ibad

          %%%%% END REPLACE_BAD_VOLS
	    end % if REPLACE_BAD_VOLS
	
	    if ROTATE_EPI
          %%%%% ROTATE_EPI
          
            % Check to see if we want a flip in z, in which case we need to
            % call fslswapdim in order for this to happen cleanly
            if sign(einfo(iexp).rotvect(9) == -1) 
              nfiles = size(flist,1);
              for ifile = 1:nfiles
                if sign(einfo(iexp).rotvect(7) == -1)  % do we want
                  % to flip x
                  xstr = '-x';
                else
                  xstr = 'x';
                end

                tmp_file = sprintf('/tmp/tmp_%s_%06d.nii', study_id, fix(rand*900000));
                % Format the FSL string
                fsl_str = sprintf('fslswapdim %s %s y -z %s', flist(ifile,:), xstr, tmp_file);
                msg = sprintf('%s\n', fsl_str);
                r = update_report(r,msg);
                status = unix(fsl_str);
                if status
                  error('Error swapping dimensions')
                end

                [fpath,fstub] = fileparts(flist(ifile,:));

                % Remove the original image file (actually, it's symlink in
                % most cases)
                unix_str = sprintf('rm %s', [fpath '/' fstub '.*']);
                msg = sprintf('%s\n', unix_str);
                r = update_report(r,msg);
                status = unix(unix_str);

                new_fname = fullfile(fpath,sprintf('%s.nii', fstub));
                % Move the tmp file in place of the original image file
                unix_str = sprintf('mv %s %s', tmp_file, new_fname);
                msg = sprintf('%s\n', unix_str);
                r = update_report(r,msg);
                status = unix(unix_str);
                if status
                  error('Error moving temporary file to run directory')
                end

                % Change the file type of the new image
                unix_str = sprintf('fslchfiletype NIFTI_PAIR %s', new_fname);
                msg = sprintf('%s\n', unix_str);
                r = update_report(r,msg);
                status = unix(unix_str);
                if status
                  error('Filetype change failed')
                end

              end
            else
                warning off  % suppress Nifti warnings
              pj_spm_reorient(cellstr(flist),einfo(iexp).rotvect);
              warning on
            end
            
          %%%%% END ROTATE_EPI
	    end % if ROTATE_EPI

        if ~exist('flist') || isempty(flist)
          if ~isempty(dir(run_outdir))
            flist = dir(run_outdir);
          elseif ~isempty(dir(run_indir))
            flist = dir(run_indir);
          end
        end
      
        if exist('flist') && ~isempty(flist)
          for ifl=1:size(flist,1)
            outdata.data{epi_idx} = ensemble_add_data_struct_row(...
                outdata.data{epi_idx},'subject_id',subid,'session',...
                isess,'ensemble_id',...
                sess.ensemble_id,'run',irun,'path',flist(ifl,:));
          end
        end

      end % if exist_epi
    end % for irun
  end % for isess
end % for isub=
