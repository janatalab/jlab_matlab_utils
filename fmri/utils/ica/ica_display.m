%ica_display
%
% batch for slice_overlay.m
%
% Tries to construct a file containing a bunch of ICA component images
%

% 8/00 Petr Janata

plot_types = {'ica_stack'};

comp2include = 155;

plot_idx = 1;
things_to_plot = plot_types(plot_idx)
nthings = length(things_to_plot);

ADD_MASK_CONTOUR = 0;

PRINT_TO_FILE = 1;
printfile = 'components.ps';
unix(['rm ' printfile]);

clear global SO
global SO

rootpath = '/data1/matlab/utils/ica/testimages/';
struct_img = '/data1/prime_long/25jan01JP/hires/n25jan01JP_004.img';


threshold_z = 0.5;			% Threshold z-score in map for
                                        % including the voxel
				
transform_types = {'sagittal','coronal','axial'};
transform_idx = 1;		% 1=sagittal, 2=coronal, 3=axial

slice_ranges = [-65 65; -60 -75; -50 70]; slice_skip = 5;
%slice_ranges = [-65 65; -60 -75; -50 70]; slice_skip = 5;

%
% Load colormaps that we might want to use
% 
%
load hot2
load hot3

newhot1 = brighten(hot,0.5);
newhot2 = brighten(hot2,0.5);
newhot3 = brighten(hot3,0.5);

my_colormap = colorcube;
my_colormap = hsv;

for ithing = 1:nthings
  disp(sprintf('Working on plots for %s', things_to_plot{ithing}))
  switch things_to_plot{ithing}
    case 'ica_stack'
      nimg = 0;
      ncomp = length(comp2include);

      SO.transform = transform_types{transform_idx};
      SO.slices = ...
	  slice_ranges(transform_idx,1):sign(diff(slice_ranges(transform_idx,:)))*slice_skip:slice_ranges(transform_idx,2);
      
      SO.labels.format = '%+02.0f';
      SO.cbar = [];
      SO.contours = [];

      SO.figure = spm_figure('GetWin', 'Graphics'); % use SPM figure window
      set(SO.figure,'DefaultLineLineWidth',1.0)

      % Get info on all the files in the directory
      P = spm_get('Files',rootpath,'*comp*.img');
  
      nvol = size(P,1);
      disp(sprintf('Found %d files in directory %s', nvol, rootpath))

      % Get the dimension info for the first file
      disp(sprintf('Getting info for desired %d files', ncomp))
      V = spm_vol(P(comp2include,:));

      % Load the structural
      nimg = nimg+1;
      SO.img(nimg).vol = spm_vol(struct_img);
      SO.img(nimg).type = 'truecolor';
      SO.img(nimg).prop = 1;
      SO.img(nimg).cmap = gray;
      
      for ic = 1:ncomp
	SO.figure = spm_figure('GetWin', 'Graphics'); % use SPM figure window
	nimg = 1;
	
	% Load the file to determine the threshold
V(ic).fname
fid = fopen(V(ic).fname,'rb','l');
	indata = fread(fid,inf,'int16')';
	fclose(fid);

	threshold = threshold_z*std(indata(:)) + mean(indata(:))
%	color_range = [threshold threshold+threshold*.5];
color_range = [1.5 1.6];     
	nimg = nimg+1;

	SO.img(nimg).vol = V(ic);
	SO.img(nimg).type = 'split';
	SO.img(nimg).range = color_range;
	SO.img(nimg).outofrange = {0, ceil(max(color_range))*2};
	SO.img(nimg).prop = 1;
	SO.img(nimg).cmap = newhot1;
	SO.cbar(end+1) = nimg;
      
	%    slice_overlay_contour  % display the slices
	slice_overlay

	if PRINT_TO_FILE
	  set(gcf,'PaperPositionMode', 'auto')
	  print('-dpsc','-painters','-noui','-append', printfile)
	end
      end % for ic=1:ncomp
  end % switch
end

