function [outdata] = ensemble_fmri_display(indata,defs)
% 
% Adapted from prime_base_display
%
% FIXME: get paths from previous steps
% 
% 02/25/06 Petr Janata
% 09/13/08 FB - adapted from autobio_display

clear global SO
global SO
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

if isfield(defs,'sinfo') && isstruct(defs.sinfo)
  sinfo=defs.sinfo;
  proc_subs = {sinfo(:).id};
  nsub_proc = length(sinfo(:));
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
check_vars = {'sinfo','plots','pathdata',{'hires','meanhires'}};
check_required_vars;

if (iscell(indata) && ~isempty(indata) && isfield(indata{1},'task') && ...
        ~isempty(strmatch('return_outdir',indata{1}.task))) || ...
        (isstruct(indata) && isfield(indata,'task') && ...
        ~isempty(strmatch('return_outdir',indata.task)))
  if exist('pathdata','var') && length(pathdata.data{1}) > 0
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

try PLOT_SINGLE_SUB = defs.display.PLOT_SINGLE_SUB;
  catch PLOT_SINGLE_SUB = 0; end
try PLOT_GROUP = defs.display.PLOT_GROUP; catch PLOT_GROUP = 0; end
try ADD_SPLIT = defs.display.ADD_SPLIT; catch ADD_SPLIT = 1; end
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
try MAX_T = defs.display.MAX_T; catch MAX_T = 10; end
try FTHRESH_MULT = defs.display.FTHRESH_MULT; catch FTHRESH_MULT = 4; end
try TTHRESH_MULT = defs.display.TTHRESH_MULT; catch TTHRESH_MULT = 3; end
try CONTOUR_INTERVAL = defs.display.COUNTOUR_INTERVAL;
  catch CONTOUR_INTERVAL = 0; end
try DEFAULT_CONTOUR_WIDTH = defs.display.DEFAULT_CONTOUR_WIDTH;
  catch DEFAULT_CONTOUR_WIDTH = 1.5; end
try USE_SPM = defs.display.USE_SPM; catch USE_SPM = 0; end
try USE_FSL = defs.display.USE_FSL; catch USE_FSL = 0; end

if USE_FSL && ~USE_SPM
  error('FSL not supported yet ...\n');
  return
elseif ~USE_FSL && ~USE_SPM
  error(['\t\tyou must specify either SPM or FSL to carry out '...
      'the analyses\n']);
  return
end

if ismember(model_id, [1:3])
  ORTHO_PMOD_MODEL = 1;
else
  ORTHO_PMOD_MODEL = 0;
end

% some path stuff
rootpath = defs.paths.outroot;

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

% Set up some general figure parameters
figdir = fullfile(rootpath,'figures',sprintf('model_%02d', model_id)); check_dir(figdir);

SO.figure = spm_figure('GetWin', 'Graphics'); % use SPM figure window
set(SO.figure,'DefaultLineLineWidth', DEFAULT_CONTOUR_WIDTH)
SO.labels.format = '%+d';
SO.labels.size = 0.1;

for iplot = 1:nplots
  plot = plots{iplot};
  disp(sprintf('Working on plots for %s', plot))
  
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

        anatdir = fullfile(sessdir, 'anatomy');
        analdir = fullfile(sessdir, 'analyses');
        spmdir = fullfile(analdir, 'spm');
        stats_dir = fullfile(spmdir, sprintf('model_%02d', model_id));

        % Clear the img stack
        SO.img = [];
        SO.cbar = [];
    	SO.contours = [];
        nimg = 0;
	
    	% Load the anatomical info
        nimg = nimg+1;
    	SO.img(nimg).type = 'truecolor';
        SO.img(nimg).prop = 1;
        SO.img(nimg).cmap = gray;

        if USE_INDIVIDUAL_HIRES
          if IS_HIRES_SPATNORM
            normprefix = 'w';
          else
            normprefix = '';
          end
          hires_fname = fullfile(anatdir,sprintf('%s%s_hires.img', normprefix,subid));
        else
          hires_fname = fullfile(spm_root,'canonical/avg152T1.nii');
        end
	
        SO.img(nimg).vol = spm_vol(hires_fname);
        SO.img(nimg).range = [];
	
        if isempty(strmatch(plot,'Hires4Subject','exact'))
	      % Load the SPM.mat file for the proper model so that we can access the
	      % xCon structure to get info about which image belongs to which contrast
	      model_fname = fullfile(stats_dir,'SPM.mat');
	    
	      msg = sprintf('Loading %s\n', model_fname);
          r = update_report(r,msg);
	      load(model_fname);

	      tdf = size(SPM.xX.X,1);  % Get the total degrees of freedom
	    
	      % Check to make sure that all the contrasts have been evaluated
	      if any(cellfun('isempty',{SPM.xCon.Vcon}))
	        msg = sprintf('Evaluating unevaluated contrasts ... \n');
            r = update_report(r,msg);
	        SPM = spm_contrasts(SPM,1:length(SPM.xCon));
	        msg = sprintf('Saving updated SPM structure: %s\n', model_fname);
            r = update_report(r,msg);
	        save(model_fname, 'SPM');
          end
	
	      xCon = SPM.xCon;
	    end % for isess
	
	    switch plot
          case 'Hires4Subject'
            ADD_MASK_CONTOUR = 0;
            % Take no further action besides plotting
          otherwise
            con_idx = strmatch(plot,{xCon.name},'exact');
	    
            if isempty(con_idx)
              msg = sprintf('No contrast matching: %s\n', plot);
              r = update_report(r,msg);
              continue
            end
	    
            threshold = get_threshold(xCon(con_idx),tdf,threshold_p_indiv);
            nimg = add2stack_split(xCon(con_idx),stats_dir,threshold);
            SO.img(nimg).cmap = newhot1;

            if ADD_CBAR, SO.cbar(end+1) = length(SO.img); end

            if ADD_SPLIT_CONTOUR
              add2stack_contour(xCon(con_idx),stats_dir,threshold,'r');
            end
        end % switch plot

        if ADD_MASK_CONTOUR
            add2stack_mask(stats_dir);
        end
	
        % Render the image
        show_it(tp)

        xfm_type = tp.transform_types{tp.transform_idx};
        % Add a title to the figure
        t = axes('position',[0 0.925 0.5 ...
	      0.05],'units','norm','visible','off');
        switch plots{iplot}
          case 'cues_primes_targets'
            title_str = sprintf(['Subject: %s\nModel: %s,\t' ...
              'F-map: Cues(red),Primes(blue),Targets(green), p-crit: %1.4f'], ...
            subid, model_proto_name, threshold_p_indiv);
          case 'motion_linear'
            title_str = sprintf(['Subject: %s\nModel: %s,\t' ...
              'F-map: Motion(red), Linear(green), p-crit: %1.4f'], ...
            subid, model_proto_name, threshold_p_indiv);
          case 'Hires4Subject'
            title_str = sprintf('%s, %s sections, %s', subid,...
              xfm_type, sinfo(isub).sessinfo(1).date);
          otherwise
            title_str = sprintf(['Subject: %s, Session %d\nModel: %s,\t%s,\tp-crit: %1.4f'], ...
              subid, isess, model_proto_name, plot, threshold_p_indiv);
        end
        text(0.1,0.5, title_str,'horizontalalign','left','fontsize',18);
        
        if PRINT_TO_FILE
          if ALLSUB_TO_ONE
            outstr = sprintf('allsubs_%s_%s',  plot, xfm_type);
            if ~strmatch(plot,'Hires4Subject','exact')
              outstr = [outstr '_' model_proto_name];
            end
            figdir = fullfile(rootpath,'figures', outstr);
            print_to_file(figdir, isub~=1);

          else
            switch plot
              case 'Hires4Subject'
                outstr = sprintf('%s_%s_%s', subid, plot, xfm_type);
              otherwise
                outstr = sprintf('%s_%s_%s', subid, model_proto_name, xfm_type);
            end

            figdir = fullfile(rootpath,'figures', outstr);
            print_to_file(figdir, ~((iplot==1)&(append_flag==0)));
            append_flag = 1;
          end
        end

        if CONVERT2PDF
          convert_to_pdf(figdir);
        end

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
    anatdir = fullfile(rootpath, 'group_analyses');
    spm_root = defs.fmri.spm.paths.spm_root;
    if exist('meanhires','var')
      hires_fname = meanhires.data{mhicol.path}{1};
    else
      hires_fname = fullfile(spm_root,'canonical/avg152T1.nii');
    end
    nimg = nimg + 1;
    SO.img(nimg).vol = spm_vol(hires_fname);
    SO.img(nimg).type = 'truecolour';
    SO.img(nimg).prop = 1;
    SO.img(nimg).cmap = gray;
    
    model_dir = fullfile(rootpath, 'analyses/spm/group', sprintf('model_%02d',model_id));
    
    stats_dir = fullfile(model_dir, plotdirstubs{iplot});
    model_fname = fullfile(stats_dir,'SPM.mat');
    
    % Load the SPM structure
    load(model_fname);
    tdf = size(SPM.xX.X,1);
    
    % Check to make sure that all the contrasts have been evaluated
    if any(cellfun('isempty',{SPM.xCon.Vcon}))
      msg = sprintf('Evaluating unevaluated contrasts ... \n');
	r = update_report(r,msg);
      SPM = spm_contrasts(SPM,1:length(SPM.xCon));
      msg = sprintf('Saving updated SPM structure: %s\n', model_fname);
	r = update_report(r,msg);
      save(model_fname, 'SPM');
    end
    
    xCon = SPM.xCon;
    
    con_idx = strmatch(plot,{xCon.name},'exact');
    threshold = get_threshold(xCon(con_idx),tdf,threshold_p_group)
	
    nimg = add2stack_split(xCon(con_idx),stats_dir,threshold);
    SO.img(nimg).cmap = newhot1;
    
    if ADD_CBAR, SO.cbar(end+1) = length(SO.img); end
    
    if ADD_SPLIT_CONTOUR
      add2stack_contour(xCon(con_idx),stats_dir,threshold,'r');
    end

    % Add the deactivations
    nimg = add2stack_split(xCon(con_idx),stats_dir,-1*threshold);
    SO.img(nimg).cmap = newhot3;

    if ADD_CBAR, SO.cbar(end+1) = length(SO.img); end
    
    if ADD_SPLIT_CONTOUR
      add2stack_contour(xCon(con_idx),stats_dir,-1*threshold,'b');
    end

    if ADD_MASK_CONTOUR
      add2stack_mask(stats_dir);
    end
    
    % Render the image
    show_it(tp)
    
    % Add a title to the figure
    t = axes('position',[0 0.925 0.5 ...
	  0.05],'units','norm','visible','off');
    group_model_str = sprintf('Group: %d subjects, Model: %s', tdf, model_proto_name);
    pstr = sprintf('p-crit: %1.4f', threshold_p_group);

    title_str = sprintf('%s\n%s\t%s', group_model_str, plot, pstr);
    
    text(0.1,0.6, title_str,'horizontalalign','left','fontsize',18);
    
    if PRINT_TO_FILE
      figdir = fullfile(rootpath,'figures', ...
	  sprintf('group_%s_%s_%s', model_proto_name,plot,tp.transform_types{tp.transform_idx}));
      print_to_file(figdir);
    end

    if PRINT_TO_FILE & ALLSUB_TO_ONE
      convert_to_pdf(figdir);
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
return  
  
function print_to_file(figdir,append)
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
return
  
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

return
  
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

return

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

return
  
function img = add2stack_contour(xCon,stats_dir,threshold,img_color)
  global SO
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
  SO.img(nimg).type = 'contour';
  SO.img(nimg).contours = ones(1,2)*threshold;
  SO.img(nimg).linespec = img_color;
  SO.contours(end+1) = nimg;

return
  
function nimg = add2stack_mask(stats_dir)
  global SO
  
  nimg = length(SO.img)+1;
  SO.img(nimg).vol = spm_vol(sprintf('%s/mask.img', stats_dir));
  SO.img(nimg).type = 'contour';
  SO.img(nimg).contours = [0.5 0.5];
  SO.img(nimg).linespec = 'w';

return
  