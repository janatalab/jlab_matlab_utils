%
% batch for slice_overlay.m
%

% 10/17/00 Petr Janata  Made compatible with Brett's 9/12/00 revision

plot_types = {'group_stats', ... 
      'indiv_stats', ...
      'group_contrast_contour' ...
      }; 

plot_idx = 1;

things_to_plot = plot_types{plot_idx}

ADD_MASK_CONTOUR = 1;
DISPLAY_ON_INDIV_BRAIN = 1;

PRINT_INDIV = 0;
PRINT_TO_FILE = 0;

clear global SO
global SO

rootpath = '/data1/ts1/';		% enter the path to where your subject
                                        % and group analysis directories lie
spm_path = '/usr/local/matlab/matlab5/toolbox/spm/spm99/';

% Load up a list of subject ids (this is used to access individual subject directories)
sinfo = get_ts1_sinfo;
subject_ids = cellstr(char(sinfo(:).id))
nsub = length(subject_ids);
proc_sub = 1:nsub;

df_indiv = 426;				% Degrees of freedom in the individual
                                        % subject models
df_group = nsub-1;

threshold_p = 0.01		% Individual subjects
                                %  2.34 = 0.01 (d.f. = 426)
                                %  3.11 = 0.001 
				%  3.75 = 0.0001
				%  4.31 = 0.00001 (d.f. = 426)

				% Group
				%  2.72 = 0.01 (d.f. = 11)
				%  4.02 = 0.001 (d.f. = 11)
				%  5.45 = 0.0001 (d.f. = 11)
				%  7.10 = 0.00001 (d.f. = 11)

max_cutoff_indiv = 10;				
			
%
% Set up slice information
%

transform_types = {'sagittal','coronal','axial'};
transform_idx = [1];		% 1=sagittal, 2=coronal, 3=axial

slice_ranges = [-65 65; -60 -75; -45 75];

slice_skip = 5;

%
% Load colormaps that we might want to use.  Note that Matthew Brett's routines
% also support several different default colormaps.
%

load hot2
load hot3

newhot1 = brighten(hot,0.5);
newhot2 = brighten(hot2,0.5);
newhot3 = brighten(hot3,0.5);

switch things_to_plot
  case 'group_stats'
    nimg = 0;
			
    threshold = abs(tinv(threshold_p, df_group))
    color_range = [threshold 8];
    
    nimg = nimg+1;
    SO.img(nimg).vol = spm_vol(sprintf('%s/group/mean.img', rootpath));
    SO.img(nimg).type = 'truecolor';
    SO.img(nimg).prop = 1;
    SO.img(nimg).cmap = gray;
    
    SO.cbar = [2 3];
    
    SO.transform = transform_types{transform_idx};
    SO.slices = ...
	slice_ranges(transform_idx,1):sign(diff(slice_ranges(transform_idx,:)))*slice_skip:slice_ranges(transform_idx,2);
    
    nanat = nimg;

    SO.labels.format = '%+d';
    SO.figure = spm_figure('GetWin', 'Graphics'); % use SPM figure window
    set(SO.figure,'DefaultLineLineWidth',2.0)

    % Set up things for the MUSIC VS REST contrast
    curdir = sprintf('%s/group/music-rest', rootpath);
    
    nimg = nanat+1;
    SO.img(nimg).vol = spm_vol(sprintf('%s/%s', curdir, 'spmT_0002.img'));
    SO.img(nimg).type = 'split';
    SO.img(nimg).range = color_range;
    SO.img(nimg).prop = 1;
    SO.img(nimg).cmap = newhot1;

    nimg = nimg+1;
    SO.img(nimg).vol = spm_vol(sprintf('%s/group/music-rest/%s', rootpath, 'spmT_0003.img'));
    SO.img(nimg).type = 'split';
    SO.img(nimg).range = color_range;
    SO.img(nimg).prop = 1;
    SO.img(nimg).cmap = newhot3;

    if ADD_MASK_CONTOUR
      nimg = nimg+1;
      SO.img(nimg).vol = spm_vol(sprintf('%s/mask.img', curdir));
      SO.img(nimg).type = 'contour';
      SO.img(nimg).cmap = hot3;
      SO.img(nimg).contours = [0.5 0.5];
      SO.img(nimg).linespec = 'g';
    end

    slice_overlay			% This is the routine that does all the work

    if PRINT_TO_FILE
      set(gcf,'PaperPositionMode', 'auto')
      print('-dpsc','-painters','-noui', sprintf('%s/group/music_vs_rest.ps',rootpath))
    end

  case 'indiv_stats'
    for isub = proc_sub
      curdir = fullfile(rootpath, subject_ids{isub},'normalized');

      threshold = abs(tinv(threshold_p, df_indiv))

      nimg = 0;
      SO.cbar = [];
      
      % Set things up for the anatomical
      
      nimg = nimg+1;
      if DISPLAY_ON_INDIV_BRAIN
	hires_dir = fullfile(rootpath, subject_ids{isub},'hires');
	eval(['ls ' hires_dir '/n*.img;']); hires_img = ans(1:end-1);
	SO.img(nimg).vol = spm_vol(hires_img);

      else % Display individual data on top of canonical image      
	SO.img(nimg).vol = spm_vol(fullfile(spm_path,'canonical/avg152T1.img'));
      end
      SO.img(nimg).type = 'truecolor';
      SO.img(nimg).prop = 1;
      SO.img(nimg).cmap = gray;
      SO.cbar(nimg) = 0;
     
      nimg = nimg+1;
      fname = sprintf('%s/spmT_0009.img', curdir)
      SO.img(nimg).vol = spm_vol(fname);
      SO.img(nimg).type = 'split';
      SO.img(nimg).cmap = hot;
      SO.img(nimg).range = [threshold max_cutoff_indiv];
      SO.cbar(nimg) = 1;
      
      SO.cbar = find(SO.cbar);

      SO.transform = transform_types{transform_idx};
      SO.slices = ...
	  slice_ranges(transform_idx,1):sign(diff(slice_ranges(transform_idx,:)))*slice_skip:slice_ranges(transform_idx,2);

      SO.figure = spm_figure('GetWin', 'Graphics'); % use SPM figure window

      slice_overlay

      if PRINT_INDIV
	set(gcf,'PaperPositionMode', 'auto')
	print('-dpsc','-painters','-noui', sprintf('%s/%s_attend.ps',rootpath,subject_ids{isub}))
      end
	
    end
    
  case 'group_contrast_contour'
    nimg = 0;
    threshold = abs(tinv(threshold_p, df_group))
    color_range = [threshold 8];
    
    nimg = nimg+1;
    SO.img(nimg).vol = spm_vol(sprintf('%s/group/mean.img', rootpath));
    SO.img(nimg).type = 'truecolor';
    SO.img(nimg).prop = 1;
    SO.img(nimg).cmap = gray;
    
    SO.transform = transform_types{transform_idx};
    SO.slices = ...
	slice_ranges(transform_idx,1):sign(diff(slice_ranges(transform_idx,:)))*slice_skip:slice_ranges(transform_idx,2);
    
    nanat = nimg;

    SO.labels.format = '%+d';
    SO.figure = spm_figure('GetWin', 'Graphics'); % use SPM figure window
    set(SO.figure,'DefaultLineLineWidth',2.0)
    
    % Set up things for the CONDITION 1 vs CONDITION 2 overlay
    nimg = nanat+1;
    SO.img(nimg).vol = spm_vol(sprintf('%s/group/t1-listen/%s', rootpath, 'spmT_0002.img'));
    SO.img(nimg).type = 'contour';
    SO.img(nimg).contours = [threshold threshold];
    SO.img(nimg).prop = 0.4;
    SO.img(nimg).cmap = newhot1;
    SO.img(nimg).linespec = 'r';
    
    nimg = nimg+1;
    SO.img(nimg).vol = spm_vol(sprintf('%s/group/t2-listen/%s', rootpath, 'spmT_0002.img'));
    SO.img(nimg).type = 'contour';
    SO.img(nimg).contours = [threshold threshold];
    SO.img(nimg).prop = 0.4;
    SO.img(nimg).cmap = newhot3;
    SO.img(nimg).linespec = 'c';

    if ADD_MASK_CONTOUR
      nimg = nimg+1;
      SO.img(nimg).vol = spm_vol(sprintf('%s/group/t1-listen/mask.img', rootpath));
      SO.img(nimg).type = 'contour';
      SO.img(nimg).cmap = hot3;
      SO.img(nimg).contours = [0.5 0.5];
      SO.img(nimg).linespec = 'g';
    end

    slice_overlay		% This is the routine that does all the work
    
    if PRINT_TO_FILE
      set(gcf,'PaperPositionMode', 'auto')
      print('-dpsc','-painters','-noui', sprintf('%s/group/t1_vs_t2_contour.ps',rootpath))
    end
end