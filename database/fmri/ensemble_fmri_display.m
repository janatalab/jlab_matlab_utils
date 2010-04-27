function [outdata] = ensemble_fmri_display(indata,defs)

% display fmri images generated within ensemble_fmri
%
%   outdata = ensemble_fmri_display(indata,defs)
% 
% REQUIRES
%   indata
%       hires
%       meanhires
%       paths
%       sinfo
%       permute_model (FIXME:have to save properly from ens_fmri_eval_perm)
%   defs
%       conj_idx
%       conjunctions
%       display
%       fmri
%       model
%       paths
%       plots
%       plotdirstubs (optional)
%       USE_PERMUTATION (requires permute_model)
%       perm_idxs (optional, requires USE_PERMUTATION) indicates indices of
%           plots to mask using permutation data
%       PERM_PROB - probability threshold for permutation mask
% 
% RETURNS
% 
% FIXME: get paths from previous steps
% FIXME: handle display using permutation testing
% 
% 02/25/06 Petr Janata
% 09/13/08 FB - adapted from autobio_display

clear global SO
global SO defaults
global CLOBBER_PS FTHRESH_MULT TTHRESH_MULT tp
global r

outdata = ensemble_init_data_struct();
outdata.type = 'fmri_display';

r = init_results_struct;

r.type = 'autobio_display';  % Identify the type of this reporting instance
r.report_on_fly = 1;

% Parse out the input data
for idata = 1:length(indata)
  if isfield(indata{idata},'type')
    switch indata{idata}.type
      case 'sinfo'
        sinfo = indata{idata};
        sinfo = sinfo.data;
        proc_subs = {sinfo(:).id};
        nsub_proc = length(proc_subs);
      case 'meanhires'
        meanhires = indata{idata};
        mhicol = set_var_col_const(meanhires.vars);
      case 'hires'
        hires = indata{idata};
        hicol = set_var_col_const(hires.vars);
      case 'paths'
        pathdata = indata{idata};
        pacol = set_var_col_const(pathdata.vars);
    end
  end
end

if isfield(defs,'plots')
  plots=defs.plots;
  nplots = length(plots);
end

% if no defs.plotdirstubs, assume plot directories are named as the plot is
if isfield(defs,'plotdirstubs')
  plotdirstubs = defs.plotdirstubs;
else
  plotdirstubs = plots;
end

if isfield(defs,'model') && isfield(defs.model,'model_id')
  model_id = defs.model.model_id;
  model_proto_name = sprintf('model_%02d', model_id);
end

% check for required vars, quit if they can't be found
check_vars = {'plots'};
check_required_vars;

% return outdir if that task is requested
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
  if ~exist('outdata','var') || (ischar(outdata) && ~exist(outdata,'dir')) ...
          || ~ischar(outdata), outdata = ''; end
  return
end

try PLOT_SINGLE_SUB = defs.display.PLOT_SINGLE_SUB;
  catch PLOT_SINGLE_SUB = 0; end
try PLOT_GROUP = defs.display.PLOT_GROUP; catch PLOT_GROUP = 0; end
try ADD_SPLIT_CONTOUR = defs.display.ADD_SPLIT_CONTOUR;
  catch ADD_SPLIT_CONTOUR = 0; end
try ADD_MASK_CONTOUR = defs.display.ADD_MASK_CONTOUR;
catch ADD_MASK_CONTOUR = 0; end
try ADD_CBAR = defs.display.ADD_CBAR; catch ADD_CBAR = 1; end
try ADD_ROI_MASKS = defs.display.ADD_ROI_MASKS; catch ADD_ROI_MASKS=0; end
try USE_INDIVIDUAL_HIRES = defs.display.USE_INDIVIDUAL_HIRES;
  catch USE_INDIVIDUAL_HIRES = 1; end
try USE_GROUP_HIRES = defs.display.USE_GROUP_HIRES;
  catch USE_GROUP_HIRES = 0; end
try IS_HIRES_SPATNORM = defs.display.IS_HIRE_SPATNORM;
  catch IS_HIRES_SPATNORM = 1; end
try PRINT_INDIV = defs.display.PRINT_INDIV; catch PRINT_INDIV = 0; end
try ALLSUB_TO_ONE = defs.display.ALLSUB_TO_ONE; catch ALLSUB_TO_ONE = 0; end
try PRINT_GROUP = defs.display.PRINT_GROUP; catch PRINT_GROUP = 1; end
try PRINT_TO_FILE = defs.display.PRINT_TO_FILE;
  catch PRINT_TO_FILE = 1; end
try CLOBBER_PS = defs.display.CLOBBER_PS; catch CLOBBER_PS = 1; end
try CONVERT2PDF = defs.display.CONVERT2PDF; catch CONVERT2PDF = 0; end
try MAX_T = defs.display.MAX_T; catch MAX_T = 0; end
try FTHRESH_MULT = defs.display.FTHRESH_MULT; catch FTHRESH_MULT = 4; end
try TTHRESH_MULT = defs.display.TTHRESH_MULT; catch TTHRESH_MULT = 3; end
try CONTOUR_INTERVAL = defs.display.COUNTOUR_INTERVAL;
  catch CONTOUR_INTERVAL = 0; end
try DEFAULT_CONTOUR_WIDTH = defs.display.DEFAULT_CONTOUR_WIDTH;
  catch DEFAULT_CONTOUR_WIDTH = 1.5; end
try USE_SPM = defs.display.USE_SPM; catch USE_SPM = 0; end
try USE_FSL = defs.display.USE_FSL; catch USE_FSL = 0; end
try ORTHO_PMOD_MODEL = defs.display.ORTHO_PMOD_MODEL;
  catch ORTHO_PMODMODEL = 0; end
try USE_PERMUTATION = defs.USE_PERMUTATION;
  catch USE_PERMUTATION = 0; end
try perm_idxs = defs.perm_idxs; catch perm_idxs = ''; end
try PERM_PROB = defs.PERM_PROB; catch PERM_PROB = 0.05; end
  
if PLOT_SINGLE_SUB
  % require subject info if plotting single subject
  check_vars = {'sinfo'};
  check_required_vars;
  if USE_INDIVIDUAL_HIRES
    if ~exist('hires','var')
      warning('no individual hires data found, using canonical!');
      USE_INDIVIDUAL_HIRES = 0;
    end
  end
end

if USE_FSL && ~USE_SPM
  error('FSL not supported yet ...\n');
elseif (~USE_FSL && ~USE_SPM) || (USE_FSL && USE_SPM)
  error(['\t\tyou must specify either SPM or FSL to carry out '...
      'the analyses\n']);
elseif USE_SPM
  if isfield(defs,'fmri') && isfield(defs.fmri,'spm') && ...
          isfield(defs.fmri.spm,'defaults')
    defaults = defs.fmri.spm.defaults;
  else
    error('USE_SPM but no defs.fmri.spm.defaults set ... \n');
  end
end

% some path stuff
rootpath = defs.paths.outroot;

if isfield(defs,'conjunctions')
  if ~isfield(defs,'conj_idx') || ~isnumeric(defs.conj_idx)
    warning('defs.conj_idx not properly specified, assuming 1\n');
    conj_idx = 1;
  else
    conj_idx = defs.conj_idx;
  end
  conjunctions = defs.conjunctions;
  CONJUNCTIONS = 1;
else
  CONJUNCTIONS = 0;
end

if isfield(defs,'display') && isstruct(defs.display) && ...
        isfield(defs.display,'tp');
  tp = defs.display.tp;
  if ~isfield(tp,'transform_types')
    tp.transform_types = {'sagittal','coronal','axial'};
  end
  if ~isfield(tp,'transform_idx')
    tp.transform_idx = 1;
  end
  xfm_type = tp.transform_types{tp.transform_idx};

  if ~isfield(tp,'non_contig_slices')
    tp.non_contig_slices = [];
  end
  if ~isfield(tp,'slice_ranges')
    tp.slice_ranges = [-70 70; 60 -90; -20 65];
  end
  if ~isfield(tp,'slice_skip')
    tp.slice_skip = 10;
  end
else
  tp.transform_types = {'sagittal','coronal','axial'};
  tp.transform_idx = 1;		% 1=sagittal, 2=coronal, 3=axial
  tp.slice_ranges = [-70 70; 60 -90; -20 65];
  tp.slice_skip = 10; % used for paper
  tp.non_contig_slices = [];
  xfm_type = tp.transform_types{tp.transform_idx};
end

try threshold_p_group = defs.display.threshold_p_group;
catch threshold_p_group = 0.01; end
try threshold_p_indiv = defs.display.threshold_p_indiv;
catch threshold_p_indiv = 0.01; end
try threshold_voxels = defs.display.threshold_voxels;
catch threshold_voxels = 0; end
try threshold_desc = defs.display.threshold_desc;
catch threshold_desc = 'none'; end

%
% Load colormaps that we might want to use
% 
%
load hot2
load hot3

newhot1 = brighten(hot,0.5); % red
newhot2 = brighten(hot2,0.5); % green/yellow
newhot3 = brighten(hot3,0.5); % blue

tmp = 0.6*hot+0.4*hot3;
newhot4 = brighten(tmp(1:35,:),0.7); % magenta

my_colormap = colorcube;
my_colormap = hsv;

conjunct = [1 0 0; 0 1 0; 1 1 0; 0 0 1; 1 0 1; 0 1 1; 1 1 1];

SO.figure = spm_figure('GetWin', 'Graphics'); % use SPM figure window
set(SO.figure,'DefaultLineLineWidth', DEFAULT_CONTOUR_WIDTH)
SO.labels.format = '%+d';
SO.labels.size = 0.1;

for iplot = 1:nplots
  plot = plots{iplot};
  disp(sprintf('Working on plots for %s', plot))
  
  if USE_PERMUTATION && (ismember(iplot,perm_idxs) || isempty(perm_idxs))
    PERM = 1;
  else
    PERM = 0;
  end

  %
  % Some plots are single subject plots, whereas others are plots of group statistics
  %
  if PLOT_SINGLE_SUB
    for isub = 1:nsub_proc      
      subid = sinfo(isub).id;
      sub_outdir = fullfile(defs.paths.inroot, subid);      
      nsess = length(sinfo(isub).sessinfo);      
      append_flag = 0;
      
      for isess = 1:nsess
        sess = sinfo(isub).sessinfo(isess);
	
        if ~sess.use_session
          msg = sprintf('\t\t\tSkipping session %d\n', isess);
          r = update_report(r,msg);
          continue
        end
	
        session_stub = sess.id;
        sessdir = fullfile(sub_outdir,session_stub);

        anatdir = fullfile(sessdir,'anatomy');
        analdir = fullfile(sessdir,'analyses');
        spmdir = fullfile(analdir,'spm');
        
        if isempty(strmatch(plot,'Hires4Subject','exact'))
          stats_dir = fullfile(spmdir,sprintf('model_%02d',model_id));
        end

        % Clear the img stack
        SO.img = [];
        SO.cbar = [];
    	SO.contours = [];
        nimg = 0;
	
    	% Load the anatomical info
        nimg = nimg+1;
    	SO.img(nimg).type = 'truecolor';
        SO.img(nimg).prop = .6;
        SO.img(nimg).cmap = gray;

        if USE_INDIVIDUAL_HIRES && exist('hires','var')
          % get individual hires
          hfilt.include.all.subject_id = {subid};
          hdata = ensemble_filter(hires,hfilt);
          hires_fname = hdata.data{hicol.path}{1};
        else
          hires_fname = fullfile(defs.fmri.spm.paths.canonical_dir,...
              'avg152T1.nii');
        end
	
        SO.img(nimg).vol = spm_vol(hires_fname);
        SO.img(nimg).range = [];
	
        % 
        %  load permutation results?
        % 
        if PERM
          % FIXME: use data from indata,'permute_models'
          %%%% HACK: locating model perm results based on expected location
          perm_fname = fullfile(stats_dir,'PermProb.hdr');
          if ~exist(perm_fname,'file')
            error('permutation results not found at %s',perm_fname);
          end
          
          Vperm = spm_vol(perm_fname);
          [Yperm,XYZperm] = spm_read_vols(Vperm);
          permmask = Yperm < PERM_PROB;
        end
        
        if ~isempty(strmatch(plot,'conjunctions','exact'))

            if ~CONJUNCTIONS || ~exist('conj_idx')
              error('conjunctions specified, but not provided in defs\n');
            end
            conj_name = conjunctions{conj_idx}.name;
            conjunction_info = conjunctions{conj_idx}.mappings;

	        model_fname = fullfile(stats_dir,'SPM.mat');
	    
	        msg = sprintf('Loading %s\n', model_fname);
            r = update_report(r,msg);
  	        load(model_fname);

            nitems = size(conjunction_info,1);
            for im = 1:nitems
              cond_name = conjunction_info{im,1};  % condition
              cont_name = conjunction_info{im,2};  % contrast name
              pval = conjunction_info{im,3};  % p-value
              color_str = conjunction_info{im,4}; % color string

              job.swd = stats_dir;
              job.title = cond_name;
              job.Im = []; % contrasts for masking
        	  job.u = pval;
          
              job.thresDesc = threshold_desc;
              job.k = threshold_voxels;
              job.Ic = strmatch(cont_name,{SPM.xCon.name},'exact');	  
      	      if isempty(job.Ic)
                msg = sprintf('No contrast matching: %s/conj\n', cont_name);
                r = update_report(r,msg);
                continue
              elseif length(job.Ic) > 1
                msg = sprintf('Many contrasts match: %s, using 1st\n', ...
                    cont_name);
                r = update_report(r,msg);
                job.Ic = job.Ic(1);
              end

              [SPM,xSPM] = spm_getSPM(job);
              
              % Initialize an output volume if this is the first volume
              if im == 1
                conj_dir = fullfile(stats_dir,'Conjunctions');
                check_dir(conj_dir);
                V.fname = fullfile(conj_dir, conj_name); 
                V.dim = SPM.xVol.DIM;
                V.mat = SPM.xVol.M;

                X = zeros([V.dim; nitems]');
              end

              % Check to see if any voxels survived the thresholding
              if isempty(xSPM.XYZ)
                continue
              end

              tmp = zeros(V.dim');
              tmp(sub2ind(SPM.xVol.DIM',xSPM.XYZ(1,:),xSPM.XYZ(2,:),xSPM.XYZ(3,:))) = 1;
              if color_str == 'w'
                X(:,:,:,im) = tmp*7;
              elseif color_str == 'k'
                X(:,:,:,im) = tmp*0;
              else
                X(:,:,:,im) = tmp*(2^(find('rgb' == color_str)-1));
              end

              nimg = length(SO.img)+1;
              vals = xSPM.Z;
              SO.img(nimg).range = [1 1.75]*xSPM.u;
              SO.img(nimg).vol = slice_overlay('blobs2vol', xSPM.XYZ, vals, xSPM.M);
              SO.img(nimg).type = 'split';
              SO.img(nimg).prop = .32;
              SO.img(nimg).contours = [0.5 1.5];

              if ADD_CBAR
                SO.cbar(end+1) = nimg;
              end

              cmap_idxs = [1:32]+8; % [1:32]
              switch color_str
                case 'r'
                  SO.img(nimg).cmap = newhot1(cmap_idxs,:);
                case 'g'
                  SO.img(nimg).cmap = newhot2(cmap_idxs,:);
                case 'b'
                  SO.img(nimg).cmap = newhot3(cmap_idxs,:);
              end

              if ADD_MASK_CONTOUR
                add2stack_mask(anatdir);
              end
            end % for im = 1:nitems

            nimg = length(SO.img)+1;
            SO.img(nimg).vol = slice_overlay('matrix2vol', sum(X,4), xSPM.M);
            SO.img(nimg).type = 'split';
            SO.img(nimg).prop = .2;
            SO.img(nimg).range = [1 2^3-1];
            SO.img(nimg).cmap = conjunct;

            if ADD_MASK_CONTOUR
              add2stack_mask(anatdir);
            end

%             if ADD_FAV_CONTOUR
%               add_fav_contour;
%             end

        elseif isempty(strmatch(plot,'Hires4Subject','exact'))
	      % Load the SPM.mat file for the proper model so that we can access the
	      % xCon structure to get info about which image belongs to which contrast
	      model_fname = fullfile(stats_dir,'SPM.mat');
	    
	      msg = sprintf('Loading %s\n', model_fname);
          r = update_report(r,msg);
	      load(model_fname);

    	  % Specify a job for spm_getSPM()
%     	  job.spmmat = {model_fname};
          job.swd = stats_dir;
  	      job.title = plot;
  	      job.Im = []; % contrasts for masking
  	      switch plot
% 	        case {'Tonreg_F'}
%               job.conspec.thresh = ...
% 	  	          finv(1-threshold_p_indiv,tonreg_eff_df,SPM.xX.erdf);
%   	        case {'PermProb'}  % Add this image directly
%     	  	  nimg = length(SO.img)+1;;
%               V = spm_vol(fullfile(stats_dir,'PermProb.img'));
%               SO.img(nimg).vol = V;
%               Y = spm_read_vols(spm_vol(V));
%               Y(Y<0.95) = NaN;
%               [Y,XYZ] = clust_filt(Y,10); % 5,10
%               SO.img(nimg).vol = slice_overlay('matrix2vol',Y,V.mat);
%               SO.img(nimg).type = 'split';
%               SO.img(nimg).range = [0.95 1];
%               SO.img(nimg).prop = 1;
%               SO.img(nimg).cmap = newhot1;
%               if ADD_CBAR, SO.cbar(end+1) = length(SO.img); end
%               continue
% 	      
%  	        case {'TonregRatio'} % Add this image directly
%                 nimg = length(SO.img)+1;
%                 if defs.fmri.perm.scale_F  % this flag is set in init_fmri_info.m
%                   fstub = 'F_Auto_NoAuto_scaled.img';
%                 else
%                   fstub = 'F_Auto_NoAuto.img';
%                 end
%                 fprintf('Loading tonreg ratio file: %s\n', fullfile(stats_dir,fstub));
%                 V = spm_vol(fullfile(stats_dir,fstub));
%                 Y = spm_read_vols(V);
%                 [Y,XYZ] = clust_filt(Y,10); % 5
%                 SO.img(nimg).vol = slice_overlay('matrix2vol',Y,V.mat);
%                 SO.img(nimg).type = 'split';
%                 SO.img(nimg).range = log10(tonreg_range);
%                 SO.img(nimg).prop = 1;
%                 SO.img(nimg).cmap = jet;
%                 if ADD_CBAR, SO.cbar(end+1) = length(SO.img); end
%                 continue
% 		
            otherwise
        	  job.u = threshold_p_indiv;
          end
          
          job.thresDesc = threshold_desc;
          job.k = threshold_voxels;
  	      job.Ic = strmatch(plot,{SPM.xCon.name},'exact');	  
  	      if isempty(job.Ic)
  	        msg = sprintf('No contrast matching: %s\n', plot);
  	        r = update_report(r,msg);
  	        continue
          elseif length(job.Ic) > 1
  	        msg = sprintf('Many contrasts match: %s, using 1st\n', plot);
  	        r = update_report(r,msg);
            job.Ic = job.Ic(1);
          end

  	      [SPM,xSPM] = spm_getSPM(job);
	  
	      % Check to see if any voxels survived the thresholding
	      if isempty(xSPM.XYZ)
	        continue
          end
	    
          if PERM
            for j=1:size(xSPM.XYZ,2)
              x=xSPM.XYZ(:,j);
              if ~permmask(x(1),x(2),x(3))
                xSPM.Z(j) = 0;
              end
            end
          end
          
	      % Add the job to the stack
	      nimg = add2stack_spm(xSPM);

	      % Add the appropriate colormap
% 	      if icont == 1
	        cmap = newhot1;
%           else
% 	        cmap = newhot3;
%           end	    
	      SO.img(nimg).cmap = cmap;

	      if ADD_CBAR, SO.cbar(end+1) = length(SO.img); end

	      if ADD_SPLIT_CONTOUR
	        add2stack_contour(xCon(con_idx),stats_dir,threshold,'r');
          end
	    end % ~if isempty(strmatch(plot,'conjunctions
	  
	    if ADD_MASK_CONTOUR
	      add2stack_mask(stats_dir);
        else
  	      ADD_MASK_CONTOUR = 0;
        end
        
	    % Render the image
	    show_it(tp)

  	    xfm_type = tp.transform_types{tp.transform_idx};
	    % Add a title to the figure
	    t = axes('position',[0 0.925 0.5 0.05],'units','norm','visible','off');
	    switch plot
	      case 'cues_primes_targets'
	        title_str = sprintf(['Subject: %s\nModel: %s,\tF-map: '...
                'Cues(red),Primes(blue),Targets(green), p-crit: %1.4f'], ...
		  	    subid, model_proto_name, threshold_p_indiv);
  	      case 'motion_linear'
  	        title_str = sprintf(['Subject: %s\nModel: %s,\t' ...
  		      'F-map: Motion(red), Linear(green), p-crit: %1.4f'], ...
  		  	  subid, model_proto_name, threshold_p_indiv);
	      case 'Hires4Subject'
	        title_str = sprintf('%s, %s sections, %s', subid,...
              xfm_type, sinfo(isub).sessinfo(1).date);
          otherwise
	        title_str = sprintf(['Subject: %s, Session %d\nModel: %s,\t%s,'...
              '\tp-crit: %1.4f'], subid, isess, model_proto_name, ...
              plot, threshold_p_indiv);
        end
	    text(0.1,0.5, title_str,'horizontalalign','left','fontsize',18);
        
	    if PRINT_TO_FILE
          if ALLSUB_TO_ONE
	        outstr = sprintf('allsubs_%s_%s',  plot, xfm_type);
	        if ~strmatch(plot,'Hires4Subject','exact')
	          outstr = [outstr '_' model_proto_name];
            end
	        figdir = fullfile(rootpath,'figures',model_proto_name,outstr);
            check_dir(fileparts(figdir));
            if PERM
              figdir = [figdir '_perm'];
            end
	        print_to_file(figdir, isub~=1);
          else
	        switch plot
	          case 'Hires4Subject'
                outstr = sprintf('%s_%s_%s',subid,plot,xfm_type);
                model_proto_name = 'hires';
              case 'conjunctions'
                outstr = sprintf('%s_conj_%d_%s_%s',model_proto_name,...
                    conj_idx,xfm_type,subid);
              otherwise
                outstr = sprintf('%s_%s_%s',subid,model_proto_name,xfm_type);
            end
	        figdir = fullfile(rootpath,'figures',model_proto_name,outstr);
            check_dir(fileparts(figdir));

            if USE_PERMUTATION
              figdir = [figdir '_perm'];
            end
            print_to_file(figdir,~((iplot==1)&(append_flag==0)));
	        append_flag = 1;
          end
          if CONVERT2PDF          
    	    convert_to_pdf(figdir);
	      end
        end % if PRINT_TO_FILE


      end % for isess=
    end % for isub
  end % if PLOT_SINGLE_SUB

  if PLOT_GROUP
    % Clear the img stack
    SO.img = [];
    SO.cbar = [];
    SO.contours = [];
    nimg = 0;
    
    % Get the anatomical for the group
    spm_root = defs.fmri.spm.paths.spm_root;
    if USE_GROUP_HIRES && exist('meanhires','var')
      hires_fname = meanhires.data{mhicol.path}{1};
    else
      hires_fname = fullfile(spm_root,'canonical/avg152T1.nii');
    end
    nimg = nimg + 1;
    SO.img(nimg).vol = spm_vol(hires_fname);
    SO.img(nimg).type = 'truecolour';
    SO.img(nimg).prop = 1;
    SO.img(nimg).cmap = gray;
    
    if isempty(strmatch(plot,'meanhires','exact'))
      model_dir = fullfile(rootpath,'analyses/spm/group',...
          sprintf('model_%02d',model_id));
    end

    if CONJUNCTIONS

        conj_name = conjunctions{conj_idx}.name;
        conjunction_info = conjunctions{conj_idx}.mappings;
        
        nitems = size(conjunction_info,1);
        for im = 1:nitems
          cond_name = conjunction_info{im,1};  % condition
          cont_name = conjunction_info{im,2};  % contrast name
          pval = conjunction_info{im,3};  % p-value
          color_str = conjunction_info{im,4}; % color string

          stats_dir = fullfile(model_dir,cond_name);
          model_fname = fullfile(stats_dir,'SPM.mat');
          
          % load the SPM structure
          load(model_fname);
          
          job.swd = stats_dir;
          job.title = cond_name;
          job.Im = []; % contrasts for masking
          job.u = pval;

          job.thresDesc = threshold_desc;
          job.k = threshold_voxels;
          job.Ic = strmatch(cont_name,{SPM.xCon.name},'exact');
          if isempty(job.Ic)
            msg = sprintf('No contrast matching: %s/conj\n', cont_name);
            r = update_report(r,msg);
            continue
          elseif length(job.Ic) > 1
            msg = sprintf('Many contrasts match: %s, using 1st\n', ...
                cont_name);
            r = update_report(r,msg);
            job.Ic = job.Ic(1);
          end

          [SPM,xSPM] = spm_getSPM(job);

          % Initialize an output volume if this is the first volume
          if im == 1
            conj_dir = fullfile(stats_dir,'Conjunctions');
            check_dir(conj_dir);
            V.fname = fullfile(conj_dir, conj_name);
            V.dim = SPM.xVol.DIM;
            V.mat = SPM.xVol.M;

            X = zeros([V.dim; nitems]');
          end

          % Check to see if any voxels survived the thresholding
          if isempty(xSPM.XYZ)
            continue
          end

          tmp = zeros(V.dim');
          tmp(sub2ind(SPM.xVol.DIM',xSPM.XYZ(1,:),xSPM.XYZ(2,:),xSPM.XYZ(3,:))) = 1;
          if color_str == 'w'
            X(:,:,:,im) = tmp*7;
          elseif color_str == 'k'
            X(:,:,:,im) = tmp*0;
          else
            X(:,:,:,im) = tmp*(2^(find('rgb' == color_str)-1));
          end

          nimg = length(SO.img)+1;
          vals = xSPM.Z;
          SO.img(nimg).range = [1 1.75]*xSPM.u;
          SO.img(nimg).vol = slice_overlay('blobs2vol', xSPM.XYZ, vals, xSPM.M);
          SO.img(nimg).type = 'split';
          SO.img(nimg).prop = 1;
          SO.img(nimg).contours = [0.5 1.5];

          if ADD_CBAR
            SO.cbar(end+1) = nimg;
          end

          cmap_idxs = [1:32]+8; % [1:32]
          switch color_str
            case 'r'
              SO.img(nimg).cmap = newhot1(cmap_idxs,:);
            case 'g'
              SO.img(nimg).cmap = newhot2(cmap_idxs,:);
            case 'b'
              SO.img(nimg).cmap = newhot3(cmap_idxs,:);
          end

          if ADD_MASK_CONTOUR
            add2stack_mask(anatdir);
          end
%           add2stack_spm(xSPM);
        end % for im = 1:nitems
        
        
    elseif isempty(strmatch(plot,'meanhires','exact'))
    
        stats_dir = fullfile(model_dir, plotdirstubs{iplot});
        model_fname = fullfile(stats_dir,'SPM.mat');

        % Load the SPM structure
        load(model_fname);

        % Specify a job for spm_getSPM()
    % 	job.spmmat = {model_fname};
        job.swd = stats_dir;
        job.title = plot;
        job.Im = []; % contrasts for masking
        job.thresDesc = threshold_desc;
        job.u = threshold_p_group;
        job.k = threshold_voxels;
        job.Ic = strmatch(plot,{SPM.xCon.name},'exact');	  
        if isempty(job.Ic)
          msg = sprintf('No contrast matching: %s\n', plot);
          r = update_report(r,msg);
          continue
        elseif length(job.Ic) > 1
          msg = sprintf('Many contrasts match: %s, using 1st\n', plot);
          r = update_report(r,msg);
          job.Ic = job.Ic(1);
        end

        [SPM,xSPM] = spm_getSPM(job);

        % Check to see if any voxels survived the thresholding
        if isempty(xSPM.XYZ)
          continue
        end

        % Add the job to the stack
        nimg = add2stack_spm(xSPM);

        % Add the appropriate colormap
    % 	if icont == 1
          cmap = newhot1;  % activation
    % 	else
    % 	  cmap = newhot3;  % deactivation
    % 	end	    
        SO.img(nimg).cmap = cmap;

        if ADD_CBAR, SO.cbar(end+1) = length(SO.img); end

        if ADD_SPLIT_CONTOUR
          add2stack_contour(xCon(con_idx),stats_dir,threshold,'r');
        end


        if ADD_MASK_CONTOUR
          add2stack_mask(group_mask_dir);  % stats_dir
        end
    
    end % if CONJUNCTIONS
    
    % Render the image
    show_it(tp)
    
    % Add a title to the figure
    t = axes('position',[0 0.925 0.5 0.05],'units','norm','visible','off');
    switch plot
      case {'PermSum'}
        group_model_str = sprintf('Sum tonreg perm, model: %s',...
            model_proto_name);
        pstr = 'Likelihood orig data by chance: 0.05';
      case {'meanhires'}
        group_model_str = sprintf('Mean Hires Image: %d subjects',...
            meanhires.data{mhicol.nsub}(1));
        pstr = '';
      otherwise
        group_model_str = sprintf('Group: %d subjects, Model: %s',...
            fix(xSPM.df(2))+1, model_proto_name);
        pstr = sprintf('p-crit: %1.4f', threshold_p_group);
    end
    
    title_str = sprintf('%s\n%s\t%s', group_model_str, plot, pstr);
    
    text(0.1,0.6, title_str,'horizontalalign','left','fontsize',18);
    
    if PRINT_TO_FILE
      if CONJUNCTIONS
        fstub = sprintf('group_%s_%s_%s',model_proto_name,conj_name,...
            tp.transform_types{tp.transform_idx});
        figdir = fullfile(conj_dir,fstub);
        check_dir(fileparts(figdir));
      elseif ~isempty(strmatch(plot,'meanhires','exact'))
        fstub = sprintf('mean_hires_%s',tp.transform_types{tp.transform_idx});
        figdir = fullfile(rootpath,'figures',fstub);
      else
        fstub = sprintf('group_%s_%s',model_proto_name,...
            tp.transform_types{tp.transform_idx});
        figdir = fullfile(rootpath,'figures',model_proto_name,fstub);
        check_dir(fileparts(figdir));
      end
      if iplot == 1, append = 0; else append = 1; end
      if PERM
        figdir = [figdir '_perm'];
      end
      print_to_file(figdir, append);
      if CONVERT2PDF
        convert_to_pdf(figdir);
      end
    end
    
  end %if PLOT_GROUP
end % for plot

r.data = SO;

%
% Various sub-functions
%

function show_it(tp)
  global SO 
  
  SO.transform = tp.transform_types{tp.transform_idx};
  
  if ~isempty(tp.non_contig_slices)
    SO.slices = tp.non_contig_slices;
  else
    SO.slices = ...
	tp.slice_ranges(tp.transform_idx,1):sign(diff(tp.slice_ranges(tp.transform_idx,:)))*tp.slice_skip:tp.slice_ranges(tp.transform_idx,2);
  end
  
  % Call routine that does the work
  slice_overlay
  
function print_to_file(figdir, append)
  global r
  if nargin < 2
    append = 0;
  end
  
  set(gcf,'PaperPositionMode', 'auto') % 
  % Make sure paper position is OK
  %	    set(gcf,'PaperUnits','inches')
  orig_pos = get(gcf,'PaperPos');
  new_pos = orig_pos - [orig_pos(1) orig_pos(2) 0 0];
  %new_pos(3:4) = new_pos(3:4)./max(new_pos(3:4));
  set(gcf,'PaperPos',new_pos);
  set(gcf,'PaperType','usletter')
  
  psname = sprintf('%s.ps', figdir);
  msg = sprintf('   Writing data to file: %s ...\n', psname);
	r = update_report(r,msg);
  if append == 0
    print('-dpsc','-painters','-noui','-adobecset', psname);
  else
    print('-append','-dpsc','-painters','-noui','-adobecset', psname);
  end

function convert_to_pdf(figdir)
  global CLOBBER_PS
  
  psname = sprintf('%s.ps', figdir);
  pdfname = sprintf('%s.pdf', figdir);
  
  output_dim = fix([850 1100]*1.0); %1.1
  unix_str = sprintf('ps2pdf -r72 -g%dx%d %s %s', output_dim(1),output_dim(2), psname, pdfname);
  unix(unix_str);
  
  % Delete the .ps file
  if CLOBBER_PS
    unix_str = sprintf('rm %s', psname);
    unix(unix_str);
  end

function threshold = get_threshold(xCon,tdf,p)
  ndf = xCon.eidf;
  switch xCon.STAT
    case 'F'
      ddf = tdf-ndf-1;
      threshold = finv(1-p,ndf,ddf);
    case 'T'
      ddf = tdf-1;
      threshold = tinv(1-p/2,ddf);
  end

function nimg = add2stack_split(xCon,stats_dir,threshold)
  global SO FTHRESH_MULT TTHRESH_MULT
  nimg = length(SO.img)+1;

  if isstr(xCon.Vspm)
    SO.img(nimg).vol = spm_vol(fullfile(stats_dir, xCon.Vspm));
  else
    if exist(xCon.Vspm.fname,'file')
      SO.img(nimg).vol = spm_vol(xCon.Vspm.fname);
    else
      SO.img(nimg).vol = spm_vol(fullfile(stats_dir, xCon.Vspm.fname));
    end
  end
  SO.img(nimg).type = 'split';
  switch xCon.STAT
    case 'F'
      SO.img(nimg).range = [threshold threshold*FTHRESH_MULT];
    case 'T'
      SO.img(nimg).range = [threshold threshold*TTHRESH_MULT];
  end
  SO.img(nimg).prop = 1;

function nimg = add2stack_contour(xCon,stats_dir,threshold,img_color, linewidth,XYZ)
  global SO
  nimg = length(SO.img)+1;

  if isstr(xCon.Vspm)
    V = spm_vol(fullfile(stats_dir, xCon.Vspm));
  else
    if exist(xCon.Vspm.fname,'file')
      V = spm_vol(xCon.Vspm.fname);
    else
      V = spm_vol(fullfile(stats_dir, xCon.Vspm.fname));
    end
  end
  
  % If we have a list of valid voxels, then we need to create a temporary
  % volume that has masked out all other voxels
  if exist('XYZ')
    Y = spm_read_vols(V);
    Yout = zeros(size(Y));
    curr_idxs = sub2ind(size(Y),XYZ(1,:)',XYZ(2,:)',XYZ(3,:)');
    Yout(curr_idxs) = Y(curr_idxs);
    Vout = V;
    Vout.fname = '/tmp/contthresh.img';
    SO.img(nimg).vol = spm_write_vol(Vout,Yout);
  else
    SO.img(nimg).vol = V;
  end
  
  SO.img(nimg).type = 'contour';
  SO.img(nimg).contours = ones(1,2)*threshold;
  SO.img(nimg).linespec = img_color;
  if exist('linewidth')
    SO.img(nimg).linewidth = linewidth;
  else
    SO.img(nimg).linewidth = 1.5;
  end
  SO.contours(end+1) = nimg;

function nimg = add2stack_mask(stats_dir)
  global SO
  
  nimg = length(SO.img)+1;
  SO.img(nimg).vol = spm_vol(sprintf('%s/mask.img', stats_dir));
  SO.img(nimg).type = 'contour';
  SO.img(nimg).contours = [0.1 0.125];
  SO.img(nimg).linespec = 'w';

function nimg = add2stack_spm(xSPM)
  global SO FTHRESH_MULT TTHRESH_MULT
  
  nimg = length(SO.img)+1;
  SO.img(nimg).vol = slice_overlay('blobs2vol', xSPM.XYZ, xSPM.Z, xSPM.M);
  SO.img(nimg).type = 'split';
  SO.img(nimg).prop = 1;

  switch xSPM.STAT
    case 'F'
      range_mult = FTHRESH_MULT;
    case 'T'
      range_mult = TTHRESH_MULT;
  end
  SO.img(nimg).range = [1 range_mult]*xSPM.u;
