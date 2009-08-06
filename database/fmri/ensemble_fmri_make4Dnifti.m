function outdata = ensemble_fmri_make4Dnifti(indata,defs)

% uses fslmerge to merge images by subject, session, or run, into a 4D nifti
%
% outdata = ensemble_fmri_make4Dnifti(indata,defs)
% 
% will merge both analyze and nifti (3d and 4d) into a new 4D nifti
% 
% file naming conventions:
%   Nsub_allruns_%dstr.nii.gz - multiple subjects
%   subid_nsess.nii.gz        - one subject, multiple sessions
%   subid_sessid_nruns.nii.gz - one subject, one session, multiple runs
% 
% output locations:
%   a subject's epi_outdir - one subject, one session, multiple runs
%   a subject's epi_outdir, first session - one subject, multiple sessions
%   exproot/group_timeseries - multiple subjects
% 
% REQUIRES
%   sinfo
%   epi data
%   path data
%   defs.make4Dnifti.CONCAT = {by_subject|by_session|by_run|all(default)}
%   defs.make4Dnifti.bet_params = string of parameters to send to bet ...
%       default is '-m -f 0.4'
%   defs.make4Dnifti.tmpdir = path to directory to place temporary files
% 
% RETURNS
%   4D nifti files in 'epi'
%   bet masks for all 4D nifti files, in 'betmasks'
% 
% NOTE: this function will z-score and brain-extract all images that it
% receives, but will not perform any other pre-processing on the images.
% Therefore, you should motion-correct, realign, slice-time-correct,
% coregister, and perform any other pre-processing on the images before
% bringing them to this script. THIS SCRIPT ASSUMES ALL IMAGES HAVE BEEN
% COREGISTERED TO A COMMON SPACE.
% 
% NOTE: even if you request to concatenate all images across subjects,
% run-level .nii, session-level .nii, and subject-level .nii files will be
% generated.
%   run-level images will contain entries in all epidata columns.
%   session-level images will contain entries in the subject_id, session,
%       and path columns
%   subject-level images will contain entries in the subject_id and path
%       columns
%   across-subject images will only contain an entry in the path column.
%   where no 'entries' are found in columns, a [] or '' will be found,
%   depending on the class of the column contents.
% 
% NOTE: only supports by_run and by_session concatenation for now
% 
% FIXME: give option to take masks from other functions, instead of
% re-calculating masks within this function
% 
% FIXME: when generating masks, make sure to follow the same convention for
% output data structs as ensemble_fmri_mask_from_{mean|norm}_epi
% 
% NOTE: assumes you're making a 4D image file from EPIs, returns paths to
% 4d images in a data struct named 'epi'
% 
% 2009.03.05 FB

global r

outdata = ensemble_init_data_struct();
outdata.type = 'make4Dnifti';

r = init_results_struct;

r.type = 'make4Dnifti';  % Identify the type of this reporting instance
r.report_on_fly = 1;

dstr = datestr(now(),30);

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
      case {'paths'}
        pathdata = indata{idata};
        pcol = set_var_col_const(pathdata.vars);
    end
  end
end

% check for required vars
check_vars = {'sinfo','pathdata','epidata'};
check_required_vars;

% return an output directory for ensemble_jobman_parallel
% for this function, it will be the path_data type specified in
% defs.output_dir_type, defaulting to epi_outdir
if (iscell(indata) && ~isempty(indata) && isfield(indata{1},'task') && ...
        ~isempty(strmatch('return_outdir',indata{1}.task))) || ...
        (isstruct(indata) && isfield(indata,'task') && ...
        ~isempty(strmatch('return_outdir',indata.task)))
  % get output directory
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
          % one output directory, save outdata = odirtype path
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
        end % if length(spatndata.data{1
      end % if ~isempty(lpathdata
    else
      if isfield(defs,'paths') && isfield(defs.paths,'grouptime') && ...
              ~isempty(defs.paths.grouptime) && exist(defs.paths.grouptime)
        outdata = defs.paths.grouptime;
      end
    end % if length(nsub_proc
  end % if exist('pathdata
  if ~exist('outdata','var') || ~exist(outdata,'dir'), outdata = ''; end
  return
end

% sinfo output data struct
outdata.vars = [outdata.vars 'sinfo'];
sinfo_idx = length(outdata.vars);
outdata.data{sinfo_idx} = ensemble_init_data_struct();
outdata.data{sinfo_idx}.type = 'sinfo';
outdata.data{sinfo_idx}.data = sinfo;

% epi output data struct
outdata.vars = [outdata.vars 'epi'];
epi_idx = length(outdata.vars);
outdata.data{epi_idx} = ensemble_init_data_struct();
outdata.data{epi_idx}.type=epidata.type;
outdata.data{epi_idx}.vars = epidata.vars;
outdata.data{epi_idx}.data{1} = {};
outdata.data{epi_idx}.data{2} = [];
outdata.data{epi_idx}.data{3} = [];
outdata.data{epi_idx}.data{4} = [];
outdata.data{epi_idx}.data{5} = {};

% mask output data struct
outdata.vars = [outdata.vars 'betmasks'];
bet_idx = length(outdata.vars);
outdata.data{bet_idx} = ensemble_init_data_struct();
outdata.data{bet_idx}.type='betmasks';
outdata.data{bet_idx}.vars = epidata.vars;
outdata.data{bet_idx}.data{1} = {};
outdata.data{bet_idx}.data{2} = [];
outdata.data{bet_idx}.data{3} = [];
outdata.data{bet_idx}.data{4} = [];
outdata.data{bet_idx}.data{5} = {};
        
% get flags
CONCAT = 'by_session';
bet_params = '-m -f 0.4 -n';
if isfield(defs,'make4Dnifti')
  if isfield(defs.make4Dnifti,'CONCAT') ...
          && ~isempty(defs.make4Dnifti.CONCAT) ...
          && ~isempty(strmatch(defs.make4Dnifti.CONCAT,...
          {'by_session','by_run'}))
    CONCAT = defs.make4Dnifti.CONCAT;
  end
  if isfield(defs.make4Dnifti,'bet_params')...
          && ischar(defs.make4Dnifti.bet_params) && ...
          ~isempty(defs.make4Dnifti.bet_params)
    bet_params = defs.make4Dnifti.bet_params;
  end
end

%%% where to store temporary files?
if isfield(defs,'make4Dnifti') && isfield(defs.make4Dnifti,'tmpdir') ...
        && ischar(defs.make4Dnifti.tmpdir) ...
        && ~isempty(defs.make4Dnifti.tmpdir)
  check_dir(defs.make4Dnifti.tmpdir);
  tmpstub = fullfile(defs.make4Dnifti.tmpdir,'tmp_make4Dnifti');
else
  tmpstub = 'tmp_make4Dnifti';
end

%%%% temporary ... must code for this
if ~isempty(strmatch(CONCAT,{'by_subject','all'}))
  error('all/by_subject concatenation not yet supported');
end

%
% START OF THE SUBJECT LOOP
%

sub_files = {};
sub_masks = {};

for isub=1:nsub_proc

  subid = sinfo(isub).id;
  msg = sprintf('\t\tPROCESSING SUBJECT (%d/%d): %s\n', isub, nsub_proc,subid);
  r = update_report(r,msg);

  % get subject epis
  sfilt.include.all.subject_id = {subid};
  sdata = ensemble_filter(epidata,sfilt);

  % get subject paths
  spdata = ensemble_filter(pathdata,sfilt);
  
  % Determine number of sessions for this subject
  nsess = length(sinfo(isub).sessinfo);

  %
  % START OF THE SESSION LOOP
  %

  sess_files = {};
  sess_masks = {};
  
  for isess = 1:nsess
    sess = sinfo(isub).sessinfo(isess);            

    if ~sess.use_session
      msg = sprintf('\t\t\tSkipping session %d\n', isess);
      r = update_report(r,msg);
      continue
    end

    % get session epis
    sessfilt.include.all.session = isess;
    sessdata = ensemble_filter(sdata,sessfilt);
    
    % get session paths
    sesspdata = ensemble_filter(spdata,sessfilt);

    %
    % START OF RUN LOOP
    %

    [runm,urun] = make_mask_mtx(sessdata.data{epicol.run});
    nruns = length(urun);

    run_files = {};
    run_masks = {};
    
    for irun = 1:nruns
      lrun = urun(irun);
      runfilt.include.all.run = lrun;
      rundata = ensemble_filter(sessdata,runfilt);

      flist = rundata.data{epicol.path};
      check_ext(flist); % sanity check ... do we have consistent extensions?
      
      runfilt.include.all.path_type = {'run_outdir'};
      pdata = ensemble_filter(pathdata,runfilt);
      outpath = pdata.data{pcol.path}{1};

      if length(flist) > 1
        % more than one file in flist, concatenate
        %%%%% WARNING: assumes all necessary within-run preprocessing
        %%%%% (realignment, coregistration, etc)
        outfstub = sprintf('%s_%d_run%d',subid,sess.ensemble_id,lrun);

        outfname = fullfile(outpath,outfstub);
        status = unix(sprintf('fslmerge -t %s %s',outfname,...
            cell2str(flist,' ')));
        if status
          error('error merging %d files into a 4D nifti',length(flist));
        end
      else
        outfname = flist{1};
      end

      % bet mask
      maskfname = fullfile(outpath,sprintf('%s_%d_run%d',subid,...
          sess.ensemble_id,lrun));
      fstr = sprintf('bet %s %s %s',outfname,maskfname,bet_params);
      status = unix(fstr);
      if status
        error('error calculating bet mask: %s, sess %d, run %d',...
            subid,isess,lrun);
      end

      run_masks{irun} = sprintf('%s_mask.nii.gz',maskfname);

      outdata.data{bet_idx} = ensemble_add_data_struct_row(...
          outdata.data{bet_idx},'subject_id',subid,'session',...
          isess,'ensemble_id',sess.ensemble_id,'run',lrun,...
          'path',maskfname);

      % calculate mean
      meanfname = fullfile(outpath,sprintf('mean_run%d',lrun));
      fstr = sprintf('fslmaths %s -Tmean %s',outfname,meanfname);
      status = unix(fstr);
      if status
        error('error calculating mean image: %s, sess %d, run %d',...
            subid,isess,lrun);
      end
      
      % calculate std
      stdfname = fullfile(outpath,sprintf('std_run%d',lrun));
      fstr = sprintf('fslmaths %s -Tstd %s',outfname,stdfname);
      status = unix(fstr);
      if status
        error('error calculating std images: %s, sess %d, run %d',...
            subid,isess,lrun);
      end
      
      % calculate z-score
      zfname = fullfile(outpath,sprintf('zscore_run%d',lrun));
      fstr = sprintf('fslmaths %s -sub %s %s',outfname,meanfname,tmpstub);
      status = unix(fstr);
      if status
        error('error mean-centering images: %s, sess %d, run %d',...
            subid,isess,lrun);
      end

      fstr = sprintf('fslmaths %s -div %s -mas %s %s',...
          tmpstub,stdfname,maskfname,zfname);
      status = unix(fstr);
      if status
        error('error calculating std images: %s, sess %d, run %d',...
            subid,isess,lrun);
      end

      unix(sprintf('rm %s*',tmpstub)); % clean up temp file

      run_files{irun} = zfname;

      outdata.data{epi_idx} = ensemble_add_data_struct_row(...
          outdata.data{epi_idx},'subject_id',subid,'session',...
          isess,'ensemble_id',sess.ensemble_id,'run',lrun,'path',zfname);

    end % for irun

    % save session-wise nifti?
    if ~isempty(strmatch(CONCAT,{'by_subject','by_session','all'}))
      % get session timeseries output path
      sessfilt.include.all.path_type = {'epi_outdir'};
      soddata = ensemble_filter(sesspdata,sessfilt);
      sodpath = soddata.data{epicol.path}{1};

      % concatenate masks
      smfname = fullfile(sodpath,sprintf('mask_concat_sess%d',isess));
      rmasks = cell2str(run_masks,' ');
      fstr = sprintf('fslmerge -t %s %s',smfname,rmasks);
      status = unix(fstr);
      if status
        error('error concatenating run masks: %s, sess %d',subid,isess);
      end
      
      fstr = sprintf('fslmaths %s -Tmean -thr 1 %s',smfname,smfname);
      status = unix(fstr);
      if status
        error('error getting mask intersect: %s, sess %d',subid,isess);
      end

      sess_masks{isess} = smfname;

      outdata.data{bet_idx} = ensemble_add_data_struct_row(...
          outdata.data{bet_idx},'subject_id',subid,'session',...
          isess,'ensemble_id',sess.ensemble_id,'run',0,'path',smfname);

      % concatenate zscored images
      szfname = fullfile(sodpath,sprintf('zscore_concat_sess%d',isess));
      rzscored = cell2str(run_files,' ');
      fstr = sprintf('fslmerge -t %s %s',szfname,rzscored);
      status = unix(fstr);
      if status
        error('error concatenating zscored runs: %s, sess %d',subid,isess);
      end
      
      % give the file a positive offset
      fstr = sprintf('fslmaths %s -add 20 -mul 200 -mas %s %s',...
          szfname,smfname,szfname);
      status = unix(fstr);
      if status
        error('error applying a positive offset to %s',szfname);
      end
      
      sess_files{isess} = szfname;
      
      outdata.data{epi_idx} = ensemble_add_data_struct_row(...
          outdata.data{epi_idx},'subject_id',subid,'session',...
          isess,'ensemble_id',sess.ensemble_id,'run',0,'path',szfname);
    end
    
  end % for isess
  
  % save subject-wise nifti
  if ~isempty(strmatch(CONCAT,{'by_subject','all'}))
    if nsess > 1
      % concatenate session masks and zscored images
    else
      sub_files{isub} = sess_files{1};
      sub_masks{isub} = sess_masks{1};
    end
  end
  
end % for isub

% save across-subject nifti
if ~isempty(strmatch(CONCAT,{'all'}))
  if nsub > 1
    % concatenate subject masks and zscored images
  else
    % no need to do anything ... files are already saved to output
    % structures, and there's nothing else to process
  end
end


% % % 
% % % subfunctions
% % % 

function xt = check_ext(flist)

% checks that extensions match between all files in a list, to avoid a
% situation where we have mixed analyze and nifti files (by extension)

nf = length(flist);
xt = '';

for i = 1:nf
  [fp,fn,fx] = fileparts(flist{i});
  if isempty(fx)
    good = false;
    error('no file extension');
  else
    if i==1
      xt = fx;
    elseif isempty(strmatch(fx,xt,'exact'))
      good = false;
      error('file extension mismatch');
    end
  end
end



% switch CONCAT
%     
%   case 'all'
%       
%     % concatenate all epis that are given
%     flist = epidata.data{epicol.path};
%     
%     % find output path
%     if nsub_proc == 1
% 
%       % one subject ... multiple sessions? multiple runs?
%       subid = proc_subs{1};
%       sfilt.include.all.subject_id = subid;
%       sfilt.include.all.path_type = 'epi_outdir';
%       ldata = ensemble_filter(pathdata,sfilt);
%       outpath = ldata.data{pcol.path}{1};
% 
%       nsess = length(sinfo.sessinfo);
%       if nsess > 1
%         outfstub = sprintf('%s_%dsess_%s',subid,nsess,dstr);
%       else
%         nruns = length(sinfo.sessinfo.use_epi_runs);
%         sessid = sinfo.sessinfo.ensemble_id;
%         outfstub = sprintf('%s_%d_%druns_%s',subid,sessid,nruns,dstr);
%       end
%     else
%       if isfield(defs.paths,'grouptime')
%         outpath = defs.paths.grouptime;
%         check_dir(outpath);
%         
%         outfstub = sprintf('%dsub_allruns_%s',nsub_proc,dstr);
%       else
%         error('no group timeseries path');
%       end
%     end
%     
%     outfname = fullfile(outpath,outfstub);
%     status = unix(sprintf('fslmerge -t %s %s',outfname,cell2str(flist,' ')));
%     if status
%       error('error merging %d files into a 4D nifti',length(flist));
%     end
%     
%     outdata.data{epi_idx} = ensemble_add_data_struct_row(...
%         outdata.data{epi_idx},'path',outfname);
%     
%   otherwise
% 
%     %
%     % START OF THE SUBJECT LOOP
%     %
% 
%     for isub=1:nsub_proc
% 
%       subid = sinfo(isub).id;
%       msg = sprintf('\t\tPROCESSING SUBJECT (%d/%d): %s\n', isub, nsub_proc,subid);
%       r = update_report(r,msg);
% 
%       sfilt.include.all.subject_id = subid;
%       sdata = ensemble_filter(epidata,sfilt);
% 
%       switch CONCAT
%         case 'by_subject'
%             
%           % concatenate all epis for a given subject
%           flist = sdata.data{epicol.path};
% 
%           sfilt.include.all.path_type = 'epi_outdir';
%           pdata = ensemble_filter(pathdata,sfilt);
%           outpath = pdata.data{pcol.path}{1};
% 
%           nsess = length(sinfo.sessinfo);
%           if nsess > 1
%             outfstub = sprintf('%s_%dsess_%s',subid,nsess,dstr);
%           else
%             nruns = length(sinfo.sessinfo.use_epi_runs);
%             sessid = sinfo.sessinfo.ensemble_id;
%             outfstub = sprintf('%s_%d_%druns_%s',subid,sessid,nruns,dstr);
%           end
% 
%           outfname = fullfile(outpath,outfstub);
%           status = unix(sprintf('fslmerge -t %s %s',outfname,...
%               cell2str(flist,' ')));
%           if status
%             error('error merging %d files into a 4D nifti',length(flist));
%           end
% 
%           if nsess > 1
%             outdata.data{epi_idx} = ensemble_add_data_struct_row(...
%                 outdata.data{epi_idx},'subject_id',subid,'path',outfname);
%           else
%             outdata.data{epi_idx} = ensemble_add_data_struct_row(...
%                 outdata.data{epi_idx},'subject_id',subid,'session',1,...
%                 'ensemble_id',sinfo.sessinfo.ensemble_id,'path',outfname);
%           end
%           
%         otherwise
% 
%           % Determine number of sessions for this subject
%           nsess = length(sinfo(isub).sessinfo);
% 
%           %
%           % START OF THE SESSION LOOP
%           %
% 
%           for isess = 1:nsess
%             sess = sinfo(isub).sessinfo(isess);            
% 
%             if ~sess.use_session
%               msg = sprintf('\t\t\tSkipping session %d\n', isess);
%               r = update_report(r,msg);
%               continue
%             end
% 
%             sessfilt.include.session = isess;
%             sessdata = ensemble_filter(sdata,sessfilt);
%             
%             switch CONCAT
%               case 'by_session'
% 
%                 % concatenate all epis for a given subject
%                 flist = sessdata.data{epicol.path};
% 
%                 sfilt.include.all.path_type = 'epi_outdir';
%                 sfilt.include.all.session = isess;
%                 pdata = ensemble_filter(pathdata,sfilt);
%                 outpath = pdata.data{pcol.path}{1};
% 
%                 outfstub = sprintf('%s_%d_%druns_%s',subid,...
%                     sess.ensemble_id,nruns,dstr);
% 
%                 outfname = fullfile(outpath,outfstub);
%                 status = unix(sprintf('fslmerge -t %s %s',outfname,...
%                     cell2str(flist,' ')));
%                 if status
%                   error('error merging %d files into a 4D nifti',length(flist));
%                 end
% 
%                 outdata.data{epi_idx} = ensemble_add_data_struct_row(...
%                     outdata.data{epi_idx},'subject_id',subid,'session',...
%                     isess,'ensemble_id',sess.ensemble_id,'path',outfname);
% 
%               otherwise
%                     
%                 %
%                 % START OF RUN LOOP
%                 %
% 
%                 [runm,urun] = make_mask_mtx(sessdata.data{epicol.run});
%                 nruns = length(urun);
% 
%                 for irun = 1:nruns
%                   lrun = urun(irun);
%                   runfilt.include.all.run = lrun;
%                   rundata = ensemble_filter(sessdata,runfilt);
%                   
%                   flist = rundata.data{epicol.path};
%                   
%                   runfilt.include.all.path_type = 'run_outdir';
%                   pdata = ensemble_filter(pathdata,runfilt);
%                   outpath = pdata.data{pcol.path}{1};
%                   
%                   outfstub = sprintf('%s_%d_run%d_%s',subid,...
%                       sess.ensemble_id,lrun,dstr);
%                   
%                   outfname = fullfile(outpath,outfstub);
%                   status = unix(sprintf('fslmerge -t %s %s',outfname,...
%                       cell2str(flist,' ')));
%                   if status
%                     error('error merging %d files into a 4D nifti',length(flist));
%                   end
% 
%                   outdata.data{epi_idx} = ensemble_add_data_struct_row(...
%                       outdata.data{epi_idx},'subject_id',subid,'session',...
%                       isess,'ensemble_id',sess.ensemble_id,'run',lrun,...
%                       'path',outfname);
%                   
%                 end % for irun
%             end % switch CONCAT
%           end % for isess
%       end % switch CONCAT
%     end % for isub=
% end % switch CONCAT
