function outdata = ensemble_fmri_eval_permute(indata,defs)

% evaluates permutations for a given fmri model
% 
%   outdata = ensemble_fmri_eval_permute(indata,defs)
% 
% after a model has been built, estimated, and had contrasts evaluated,
% this function will evaluate permutations to establish a null distribution
% based on the model data.
% 
% NOTE: currently only supports permutation of SPM models
% 
% REQUIRES
%   indata
%       sinfo
%       modelspec - paths to model spec files (SPM.mat, or design.fsf)
% 
% OPTIONAL INPUTS
%   indata
%       modelspec - contains a pre-built model, if you want to build in one
%           step and then estimate in another step
%       epidata - the epi data used to create the model(s) in modelspec.
%           These are not directly used by this function, but if you are
%           using copy_to_local, including this indata will ensure that the
%           proper EPI files for permutation estimation are also copied to
%           local disk with the model file.
%   defs.fmri
%   defs.model
%   defs.model.permute_params
% 
% OUTPUT
%   sinfo
%   modelspec
%   permutespec
% 
% FIXME: no permutation testing handled in FSL right now ....
% 
% FB 2008.08.27

global r defaults

spm_defaults

outdata = ensemble_init_data_struct();
outdata.type = 'permute_model';

r = init_results_struct;

r.type = 'fmri_permute_model';  % Identify the type of this reporting instance
r.report_on_fly = 1;

% Parse out the input data
for idata = 1:length(indata)
  if isfield(indata{idata},'type')
      switch indata{idata}.type
%         case 'sinfo'
%           sinfo = indata{idata}.data;
%           proc_subs = {sinfo(:).id};
        case 'paths'
          pathdata = indata{idata};
          pacol = set_var_col_const(pathdata.vars);
        case 'modelspec'
          modelspec = indata{idata};
          mocol = set_var_col_const(modelspec.vars);
      end
  end
end

% check for required vars, quit if they can't be found
check_vars = {'modelspec','pathdata'};
check_required_vars;

% nsub_proc = length(proc_subs);

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
        sfilt.include.all.path_type = {'anal_outdir'};
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

nmod = length(modelspec.data{mocol.path});
for imod=1:nmod

    % get model
    mid = modelspec.data{mocol.model_id}(imod);
    midx = find(mid == [defs.model.model_id]);
    try curr_model = defs.model(midx);
      model_id = curr_model.model_id;
    catch
      msg = sprintf('error finding model for model id %d,SKIPPING\n',mid);
      r = update_report(r,msg);
      continue
    end

    % get model specification
    mfile = modelspec.data{mocol.path}{imod};
    if ~exist(mfile,'file')
      msg = sprintf('can not find model file %s, SKIPPING\n',mfile);
      r = update_report(r,msg);
      continue
    end

    % get paths for this model's subject
    subid = modelspec.data{mocol.subject_id}{imod};
    pf.include.all.path_type = {'sub_indir'};
    pf.include.all.subject_id = {subid};
    pf.include.all.session = modelspec.data{mocol.session}(imod);
    pdata = ensemble_filter(pathdata,pf);
    if ~isempty(pdata.data{pacol.path})
      figdir = fullfile(pdata.data{pacol.path}{1});
    else
      figdir = fullfile(defs.paths.outroot,'figures');
    end
    check_dir(figdir);
    
    % evaluate permutations
    %%%%% ASSUMES SPM
    if isfield(curr_model,'permute_params')

      % get permutation parameters
      pp = curr_model.permute_params;
      try maxiter=pp.maxiter;
      catch
        if isfield(defs.model,'maxiter')
          maxiter=defs.model.maxiter;
        else
          maxiter=50;
        end
      end
      
      est_method = curr_model.est_method;
      switch est_method
        case {'Classical'}
          pjobs{1}.stats{1}.fmri_est.method.Classical = 1;
        case {'Bayesian'}
          pjobs{1}.stats{1}.fmri_est.method.Bayesian = curr_model.bayes;
        case {'Bayesian2'}
          pjobs{1}.stats{1}.fmri_est.method.Bayes2 = 1;
        otherwise
          msg = sprintf('ERROR: Unknown estimation method: %s\n', est_method);
          r = update_report(r,msg);
          continue
      end % switch est_method
      if isfield(pp,'stopping_rule_fun')
        stoprulefh = parse_fh(pp.stopping_rule_fun);
      else
        stoprulefh = '';
      end
      msg = sprintf('executing a maximum of %d iterations',maxiter);
      r = update_report(r,msg);

      % get design matrix figure file handle
      pp.figpath = fullfile(figdir,sprintf('%s_permute_desmat_model%d.ps',...
          subid,mid));
      if exist(pp.figpath,'file'), delete(pp.figpath), end
      
      % get SPM struct from job
      origSPM = load(mfile);
      origSPM = origSPM.SPM;
      tempSPM = rmfield(origSPM,{'xVol','Vbeta','VResMS','VM','xCon','swd'});
      model_outdir = fileparts(mfile);

      nvox = size(origSPM.xVol.XYZ,2);
      pp.XYZ = origSPM.xVol.XYZ;
      pp.vol_dim = origSPM.xVol.DIM';

      % load residuals from original job
      srcres_fname = fullfile(model_outdir,'ResMS.img');
      if ~exist(srcres_fname,'file')
        msg = sprintf('residuals not found for model %d, SKIPPING\n',mid);
        r = update_report(r,msg);
        continue
      end
      Vsrc = spm_vol(srcres_fname);  % map the filename
      fprintf(1,'Reading source data: %s\n', srcres_fname);
      Ysrc = spm_read_vols(Vsrc); % read the data

      % create empty permute residual struct
      Yperm = zeros([size(Ysrc) 0]);

      % init a place for permutation history
      pp.history.std_mu = [];
      pp.history.std_sig = [];
      pp.history.std_prob = [];
      
      % iterate until stopping rule or max iter
      done = false;
      iter = 0;
      permutations = []; % FIXME: save perms to make sure none are repeated ??
      figure();
      while (~done)
        iter = iter + 1;
        curr_dir = fullfile(model_outdir,sprintf('perm%04d',iter));
        check_dir(curr_dir);
        delete(fullfile(curr_dir,'*'));

        % permute, save data
        [SPM,pp] = fmri_permute_desmat(tempSPM,pp);
        outfname = fullfile(curr_dir,'SPM.mat');
        save(outfname,'SPM');
        
        % estimate model
        pjobs{1}.stats{1}.fmri_est.spmmat = {outfname};
        warning off
        msg = sprintf('Launching permutation %d ...\n',iter);
        r = update_report(r,msg);
        spm_jobman('run',pjobs);
        warning on

        % load new residuals
        permres_fname = fullfile(curr_dir,'ResMS.img');
        Vperm(iter) = spm_vol(permres_fname);
        Yperm(:,:,:,iter) = spm_read_vols(Vperm(iter));

        % estimate null distribution
        [Ymu,Ysig,Yprob,YpermZ] = fmri_permute_eval_dist(Ysrc,Yperm,pp.XYZ);

        % evaluate stopping rule
        pp.Ymu = Ymu;
        pp.Ysig = Ysig;
        pp.Yprob = Yprob;
        pp.YpermZ = YpermZ;
        if ~isempty(stoprulefh)
          [done,pp] = stoprulefh(pp);
        end

        % save this permutation's parameters
        pp.lastYmu = pp.Ymu;
        pp.lastYsig = pp.Ysig;
        pp.lastYprob = pp.Yprob;
        pp.history.std_mu(iter) = nanstd(Ymu(:));
        pp.history.std_sig(iter) = nanstd(Ysig(:));
        pp.history.std_prob(iter) = nanstd(Yprob(:));

        if iter == 1, continue, end

        % plot mean and standard deviation of this perm's params
        subplot(3,1,1);
        plot(pp.history.std_mu(2:end));
        title(sprintf('std of Ymu over %d iterations',iter));
        subplot(3,1,2);
        plot(pp.history.std_sig(2:end));
        title(sprintf('std of Ysig over %d iterations',iter));
        subplot(3,1,3);
        plot(pp.history.std_prob(:,2:end));
        title(sprintf('std of Yprob over %d iterations',iter));

        if iter == maxiter, done = true; end
      end

      % Write the probability volume to file
      Vprob = Vsrc;
      Vprob.fname = fullfile(model_outdir,'PermProb.img');
      Vprob.descrip = 'Probability that ResMS is lower than permuted models';
%       if isfield(Vprob,'private'), Vprob = rmfield(Vprob,'private'); end

      Vprob = spm_create_vol(Vprob);
      Vprob = spm_write_vol(Vprob,Yprob);  

      %%%% FIXME: make permutation report, including model stats and status
      
      % print std plots of permutation parameters
      figoutpath = fullfile(figdir,sprintf('%s_permute_progress_model%d.ps',...
          subid,mid));
      print('-dpsc','-painters','-noui','-adobecset',figoutpath);

    else
        msg = sprintf('no model permutation params provided for model %d,SKIPPING\n',...
            modelspec.data{mocol.model_id}(imod));
        r = update_report(r,msg);
        continue
    end % if isfield(curr_model,'permute_params
end % for imod=1:nmod
