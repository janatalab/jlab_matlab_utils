function varargout = slice_overlay(action, varargin);
% Function to display + manage slice display 
% Slice display works on a global structure SO
% with fields 
%  - img - array of images to display
%        - img structs contain fields
%             vol - vol struct info (see spm_vol)
%             cmap - colormap for this image
%             nancol - color for NaN. If scalar, this is an index into
%                    the image cmap.  If 1x3 vector, it's a colour
%             prop - proportion of intensity for this cmap/img 
%                    if = Inf, gives split cmap effect where values of 
%                    this cmap override previous image cmap values
%             func - function to apply to image before scaling to cmap
%                    (and therefore before min/max thresholding. E.g. a func of
%                    'i1(i1==0)=NaN' would convert zeros to NaNs
%             range - 2x1 vector of values for image to distribute colormap across
%                    the first row of the colormap applies to the first
%                    value in 'range', and the last value to the second
%                    value in 'range'
%             outofrange - behavior for image values to the left and
%                    right of image limits in 'range'.  Left means
%                    colormap values < 1, i.e for image values <
%                    range(1), if (range(1)<range(2)), and image values >
%                    range(1) where (range(1)>range(2)). If missing,
%                    display min (for Left) and max (for Right) value from colormap. 
%                    Otherwise should be a 2 element cell array, where
%                    the first element is the colour value for image values
%                    left of 'range', and the second is for image values
%                    right of 'range'.  Scalar values for
%                    colour index the colormap, 3x1 vectors are colour
%                    values.  An empty array attracts default settings
%                    appropriate to the mode - i.e. transparent colour (where 
%                    SO.prop ~= Inf), or split colour.  Empty cells
%                    default to 0. 0 specifies that voxels with this
%                    colour do not influence the image (split =
%                    background, true = black)
%            hold  - resampling order for image (see spm_sample_vol) -
%                    default 1
%            background - value when resampling outside image - default
%                    NaN
%            
% - transform - either - 4x4 transformation to apply to image slice position,
%             relative to mm given by slicedef, before display
%               or     - text string, one of axial, coronal, sagittal
%                        These orientations assume the image is currently
%                        (after its mat file has been applied) axially
%                        oriented
% - slicedef - 2x3 array specifying dimensions for slice images in mm
%             where rows are x,and y of slice image, and cols are neg max dim,
%             slice separation and pos max dim
% - slices   - vector of slice positions in mm in z (of transformed image)
% - figure    - figure handle for slice display figure
% - area      struct with fields
%                  position - bottom left, x size y size 1x4 vector of
%                      area in which to display slices
%                  units    - one of
%                    inches,centimeters,normalized,points,{pixels}
%                  halign - one of left,{center},right
%                  valign - one of top,{middle},bottom
% - xsliceno  - no of slices to display across figure (defaults to an optimum)
% - refreshf  - flag - if set or empty, refresh axis info for figure
%             else assume this is OK
% - cbar      - if empty, missing, no colourbar.  If an array of integers, then 
%             indexes img array, and makes colourbar for each cmap for
%             that img.  Cbars specified in order of appearance L->R
% - labels - struct can be absent (-> default numerical labels)
%                  empty (SO.labels = []) (no labels) or contain fields 
%                  colour - colour for label text 
%                  size - font size in units normalized to slice axes 
%                  format - if = cell array of strings =
%                  labels for each slice in Z.  If is string, specifies
%                  sprintf format string for labelling in distance of the
%                  origin (Xmm=0, Ymm=0) of each slice from plane containing
%                  the AC, in mm, in the space of the transformed image
%  
%  
% V 0.7 1/8/00  
% More or less  beta - take care.  Please report problems to  
% Matthew Brett matthew@mrc-cbu.cam.ac.uk

global SO

if nargin < 1
  action = 'display';
else
  action = lower(action);
end

switch action
 case 'checkso'
  checkso;
 case 'getcmap'
  varargout = {getcmap(varargin{1})};
 case 'volmaxmin'
  [mx mn] = volmaxmin(varargin{1});
  varargout = {mx, mn};
 case 'display'

% check and/or initialise SO struct
checkso;

% get coordinates for plane
X=1;Y=2;Z=3;
dims = SO.slicedef;
xmm = dims(X,1):dims(X,2):dims(X,3);
ymm = dims(Y,1):dims(Y,2):dims(Y,3);
zmm = SO.slices;
[y x] = meshgrid(ymm,xmm');
vdims = [length(xmm),length(ymm),length(zmm)];

% no of slices, and panels (an extra for colorbars)
nslices = vdims(Z);
minnpanels = nslices;
cbars = 0;
if is_there(SO,'cbar')
  cbars = length(SO.cbar);
  minnpanels = minnpanels+cbars;
end

% get figure data
% if written to, the axes may be specified already
figno = figure(SO.figure);

% (re)initialize axes and stuff

% check if the figure is set up correctly
if ~SO.refreshf
  axisd = flipud(findobj(SO.figure, 'Type','axes','Tag', 'slice overlay panel'));
  npanels = length(axisd);
  if npanels < vdims(Z)+cbars;
    SO.refreshf = 1;
  end
end
if SO.refreshf
  % clear figure and axis store
  clf
  axisd = [];

  % prevent print inversion problems
  set(figno,'InvertHardCopy','off');
  
  % calculate area of display in pixels
  parea = SO.area.position;
  if ~strcmp(SO.area.units, 'pixels')
    ubu = get(SO.figure, 'units');
    set(SO.figure, 'units','pixels');
    tmp = get(SO.figure, 'Position');
    ascf = tmp(3:4);
    if ~strcmp(SO.area.units, 'normalized')
      set(SO.figure, 'units',SO.area.units);
      tmp = get(SO.figure, 'Position');
      ascf = ascf ./ tmp(3:4);
    end
    set(figno, 'Units', ubu);
    parea = parea .* repmat(ascf, 1, 2);
  end
  asz = parea(3:4);
  
  % by default, make most parsimonious fit to figure
  yxratio = length(ymm)*dims(Y,2)/(length(xmm)*dims(X,2));
  if ~is_there(SO, 'xslices')
    % iteration needed to optimize, surprisingly.  Thanks to Ian NS
    axlen(X,:)=asz(1):-1:1;
    axlen(Y,:)=yxratio*axlen(X,:);
    panels = floor(asz'*ones(1,size(axlen,2))./axlen);
    estnpanels = prod(panels);
    tmp = find(estnpanels >= minnpanels);
    if isempty(tmp)
      error('Whoops, cannot fit panels onto figure');
    end
    b = tmp(1); % best fitting scaling
    panels = panels(:,b);
    axlen = axlen(:, b);
  else
    % if xslices is specified, assume X is flush with X figure dimensions
    panels(X) = SO.xsliceno;
    axlen(X) = asz(X)/panels(X);
  end
  
  % Axis dimensions are in pixels.  This prevents aspect ratio rescaling
  panels(Y) = ceil(minnpanels/panels(X));
  axlen(Y) = axlen(X)*yxratio;
  
  % centre (etc) panels in display area as required
  divs = [Inf 2 1];the_ds = [0;0];
  the_ds(X) = divs(strcmp(SO.area.halign, {'left','center','right'}));
  the_ds(Y) = divs(strcmp(SO.area.valign, {'bottom','middle','top'}));
  startc = parea(1:2)' + (asz'-(axlen.*panels))./the_ds;
  
  % make axes for panels
  r=0;c=1;
  npanels = prod(panels);
  lastempty = npanels-cbars;
  for i = 1:npanels
    % panel userdata
    if i<=nslices
      u.type = 'slice';
      u.no   = zmm(i);
    elseif i > lastempty
      u.type = 'cbar';
      u.no   = i - lastempty;
    else
      u.type = 'empty';
      u.no   = i - nslices;
    end
    axpos = [r*axlen(X)+startc(X) (panels(Y)-c)*axlen(Y)+startc(Y) axlen'];
    axisd(i) = axes(...
	'Parent',figno,...
	'XTick',[],...
	'XTickLabel',[],...
	'YTick',[],...
	'YTickLabel',[],...
	'Box','on',...
	'XLim',[1 vdims(X)],...
	'YLim',[1 vdims(Y)],...
	'Units', 'pixels',...
	'Position',axpos,...
	'Tag','slice overlay panel',...
	'UserData',u);
    r = r+1;
    if r >= panels(X)
      r = 0;
      c = c+1;
    end
  end
end

% sort out labels
if is_there(SO,'labels')
  labels = SO.labels;
  if iscell(labels.format)
    if length(labels.format)~=vdims(Z)
      error(...
	  sprintf('Oh dear, expecting %d labels, but found %d',...
		  vdims(Z), length(labels.contents)));
    end
  else
    % format string for mm from AC labelling
    fstr = labels.format;
    labels.format = cell(vdims(Z),1);
    acpt = SO.transform * [0 0 0 1]';
    for i = 1:vdims(Z)
      labels.format(i) = {sprintf(fstr,zmm(i)-acpt(Z))};
    end
  end
end

% modify colormaps with any new colours
nimgs = length(SO.img);
lrn = zeros(nimgs,3);
cmaps = cell(nimgs);
for i = 1:nimgs
  cmaps(i)={SO.img(i).cmap};
  lrnv = {SO.img(i).outofrange{:}, SO.img(i).nancol};
  for j = 1:length(lrnv)
    if prod(size(lrnv{j}))==1
      lrn(i,j) = lrnv{j};
    else
      cmaps(i) = {[cmaps{i}; lrnv{j}(1:3)]};
      lrn(i,j) = size(cmaps{i},1);
    end
  end
end

% cycle through slices displaying images
nvox = prod(vdims(1:2));
pandims = [vdims([2 1]) 3]; % NB XY transpose for display

zimg = zeros(pandims);
for i = 1:nslices
  ixyzmm = [x(:)';y(:)';ones(1,nvox)*zmm(i);ones(1,nvox)];
  img = zimg;
  
  % construct overlay
  for j = 1:nimgs    
    if ~any(ismember(j,SO.contours))
      thisimg = SO.img(j);
      
      % to voxel space of image
      vixyz = inv(SO.transform*thisimg.vol.mat)*ixyzmm;
      % raw data 
      i1 = spm_sample_vol(thisimg.vol,vixyz(X,:),vixyz(Y,:),vixyz(Z,:), ...
	  [thisimg.hold thisimg.background]);
      if is_there(thisimg, 'func')
	eval(thisimg.func);
      end
      % transpose to reverse X and Y for figure
      i1 = reshape(i1, vdims(1:2))';
      % rescale to colormap
      [csdata badvals]= scaletocmap(...
	  i1,...
	  thisimg.range(1),...
	  thisimg.range(2),...
	  cmaps{j},...
	  lrn(j,:));
      % take indices from colormap to make true colour image
      iimg = reshape(cmaps{j}(csdata(:),:),pandims);
      tmp = repmat(logical(~badvals),[1 1 3]);
      if thisimg.prop ~= Inf % truecolor overlay
	img(tmp) = img(tmp) + iimg(tmp)*thisimg.prop;
      else % split colormap effect
	img(tmp) = iimg(tmp);
      end
    end
  end

  % threshold out of range values
  img(img>1) = 1;
  
  image('Parent', axisd(i),...
	'CData',img);
  if is_there(SO,'labels')
    text('Parent',axisd(i),...
	'Color', labels.colour,...
	'FontUnits', 'normalized',...
	'VerticalAlignment','bottom',...
	'HorizontalAlignment','left',...
	'Position', [1 1],...
	'FontSize',labels.size,...
	'String', labels.format{i});
  end
   
  % add contours if desired
  for j = 1:nimgs 
    if any(ismember(j,SO.contours))
      thisimg = SO.img(j);
      % to voxel space of image
      vixyz = inv(SO.transform*thisimg.vol.mat)*ixyzmm;
      % raw data 
      i1 = spm_sample_vol(thisimg.vol,vixyz(X,:),vixyz(Y,:),vixyz(Z,:), ...
	  [thisimg.hold thisimg.background]);
      if is_there(thisimg, 'func')
	eval(thisimg.func);
      end
      % transpose to reverse X and Y for figure
      i1 = reshape(i1, vdims(1:2))';
      
      % add contour plot
      axes(axisd(i));
      set(axisd(i),'NextPlot','add');
      contour(i1,[thisimg.range(1) thisimg.range(1)],thisimg.linestr)
    end
  end
end

for i = (nslices+1):npanels
   set(axisd(i),'Color',[0 0 0]);
end
% add colorbar(s) 
for i = 1:cbars
  axno = axisd(end-cbars+i);
  cbari = SO.img(SO.cbar(i));
  cml = size(cbari.cmap,1);
  p = get(axno, 'Position');; % position of last axis
  cw = p(3)*0.2;
  ch = p(4)*0.75;
  pc = p(3:4)/2;
  [axlims idxs] = sort(cbari.range);
  a=axes(...
      'Parent',figno,...
      'XTick',[],...
      'XTickLabel',[],...
      'Units', 'pixels',...
      'YLim', axlims,...   
      'FontUnits', 'normalized',...
      'FontSize', 0.075,...
      'YColor',[1 1 1],...
      'Tag', 'cbar',...
      'Box', 'off',...
      'Position',[p(1)+pc(1)-cw/2,p(2)+pc(2)-ch/2,cw,ch]...
      );
  ih = image('Parent', a,...
	'YData', axlims(idxs),...     
	'CData', reshape(cbari.cmap,[cml,1,3]));

end

% end switch action
end

return

function checkso
% checks and fills SO structure
global SO

% figure
if is_there(SO, 'figure')
  try
    if ~strcmp(get(SO.figure,'Type'),'figure')
      error('Figure handle is not a figure')
    end
  catch
    error('Figure handle is not a valid figure')
  end
else
  % no figure handle. Try spm figure, then gcf
  SO.figure = spm_figure('FindWin', 'Graphics'); 
  if isempty(SO.figure)
    SO.figure = gcf;
  end
end
% set defaults for SPM figure 
if strcmp(get(SO.figure, 'Tag'),'Graphics')
  SO.area.position = [0 0 1 0.92]; % position figure nicely for SPM
  SO.area.units = 'normalized';
  SO.area.valign = 'top';
end

% orientation; string or 4x4 matrix
orientn = [];
SO = set_def(SO, 'transform', 'axial');
if ischar(SO.transform)
  orientn = find(strcmpi(SO.transform, {'axial','coronal','sagittal'}));
  if isempty(orientn)
    error(sprintf('Unexpected orientation %s', SO.transform));
  end
  ts = [0 0 0 0 0 0 1 1 1;...
      0 0 0 pi/2 0 0 1 -1 1;...
      0 0 0 pi/2 0 -pi/2 -1 1 1];
  SO.transform = spm_matrix(ts(orientn,:));
end
% default slice size, slice matrix depends on orientation
if ~is_there(SO,'slicedef' | ~is_there(SO, 'slices'))
  % take image sizes from first image
  V = SO.img(1).vol;
  D = V.dim(1:3);
  T = SO.transform * V.mat;
  vcorners = [1 1 1; D(1) 1 1; 1 D(2) 1; D(1:2) 1; ...
	     1 1 D(3); D(1) 1 D(3); 1 D(2:3) ; D(1:3)]';
  corners = T * [vcorners; ones(1,8)];
  SC = sort(corners');
  vxsz = sqrt(sum(T(1:3,1:3).^2));
  
  SO = set_def(SO, 'slicedef',...
    [SC(1,1) vxsz(1) SC(8,1);SC(1,2) vxsz(2) SC(8,2)]);
  SO = set_def(SO, 'slices',[SC(1,3):vxsz(3):SC(8,3)]);
end

% no colourbars by default
SO = set_def(SO, 'cbars', []);

% always refresh figure window, by default
SO = set_def(SO, 'refreshf', 1);  

% labels
defstruct = struct('colour',[1 1 1],'size',0.075,'format', '%+30f');
if ~isfield(SO, 'labels') % no field, -> default
  SO.labels = defstruct;
elseif ~isempty(SO.labels) % empty -> no labels
  % colour for slice labels
  SO.labels = set_def(SO.labels, 'colour', defstruct.colour); 
  % font size normalized to image axis
  SO.labels = set_def(SO.labels, 'size', defstruct.size); 
  % format string for slice labels
  SO.labels = set_def(SO.labels, 'format', defstruct.format); 
end
					 
% figure area stuff
defarea = struct('position',[0 0 1 1],'units','normalized');
SO = set_def(SO, 'area', defarea);
if ~is_there(SO.area, 'position')
  SO.area = defarea;
end
if ~is_there(SO.area,'units')
  if (all(SO.area.position>=0 & SO.area.position<=1))
    SO.area.units = 'normalized';
  else
    SO.area.units = 'pixels';
  end
end
SO.area = set_def(SO.area,'halign', 'center');
SO.area = set_def(SO.area,'valign', 'middle');

% fill various img arguments
% would be nice to use set_def, but we can't

% set colour intensities as we go
remcol = 1;
for i = 1:length(SO.img)
  if ~is_there(SO.img(i),'hold')
    SO.img(i).hold = 1;
  end
  if ~is_there(SO.img(i),'background')
    SO.img(i).background = NaN;
  end
  if ~is_there(SO.img(i),'prop')
    % default is true colour
    SO.img(i).prop = remcol/(length(SO.img)-i+1);
    remcol = remcol - SO.img(i).prop;
  end
  if ~is_there(SO.img(i),'range')
    [mx mn] = volmaxmin(SO.img(i).vol);
    SO.img(i).range = [mn mx];
  end
  if ~is_there(SO.img(i),'cmap')
    if SO.img(i).prop == Inf; % split map
      SO.img(i).cmap = check_map(i, 'spm_cols');
    else                  % true colour
      SO.img(i).cmap = check_map(i, 'actc');
    end
  end  
  if ~is_there(SO.img(i),'outofrange')
    % this can be complex, and depends on split/true colour
    if SO.img(i).prop == Inf % split colour
      if xor(SO.img(i).range(1) < SO.img(i).range(2), ...
	     SO.img(i).range(2) < 0)
	SO.img(i).outofrange = {[0],size(SO.img(i).cmap,1)};
      else
	SO.img(imgno).outofrange={[1], [0]};
      end
    else            % true colour
      SO.img(i).outofrange = {1,size(SO.img(i).cmap,1)};
    end
  end
  for j=1:2
    if isempty(SO.img(i).outofrange{j})
      SO.img(i).outofrange(j) = {0};
    end
  end
  if ~is_there(SO.img(i),'nancol')
    SO.img(i).nancol = 0;
  end
end  
return

function tf = is_there(a, fname)
% returns true if field fname is present in struct a, and not empty
tf = isfield(a, fname);
if tf
  tf = ~isempty(getfield(a, fname));
end
return

function [img, badvals]=scaletocmap(inpimg,mn,mx,cmap,lrn)
img = (inpimg-mn)/(mx-mn);  % img normalized to mn=0,mx=1
cml = size(cmap,1);
if cml==1 % values between 0 and 1 -> 1
  img(img>=0 & img<=1)=1;
else
  img = img*(cml-1)+1;
end
outvals = {img<1, img>cml, isnan(img)};
img= round(img);
badvals = zeros(size(img));
for i = 1:length(lrn)
  if lrn(i)
    img(outvals{i}) = lrn(i);
  else
    badvals = badvals | outvals{i};
    img(outvals{i}) = 1;
  end    
end
return

function st = set_def(st, fld, def)
if ~is_there(st, fld)
  st = setfield(st, fld, def);
end

function cm = check_map(idx, mapstr)
try
  load(mapstr)
catch
  error(sprintf(['No colormap defined for image %d, '...
		 'and default map %s is not on the path'],...
		idx, mapstr));
end
cm = eval(mapstr);
return

function [mx,mn] = volmaxmin(vol)
mx = -Inf;mn=Inf;
for i=1:vol.dim(3),
  tmp = spm_slice_vol(vol,spm_matrix([0 0 i]),vol.dim(1:2),[0 NaN]);
  tmp = tmp(find(finite(tmp(:))));
  if ~isempty(tmp)
    mx = max([mx; tmp]);
    mn = min([mn; tmp]);
  end
end
return

function cmap = getcmap(acmapname)
% get colormap of name acmapname
if ~isempty(acmapname)
  cmap = evalin('base',acmapname,'[]');
  if isempty(cmap) % not a matrix, is it...
    % a colour name?
    tmp = strcmp(acmapname, {'red','green','blue'});
    if any(tmp)
      cmap = zeros(64,3);
      cmap(:,tmp) = ((0:63)/63)';
    else
      % a .mat file?
      [p f e] = fileparts(acmapname);
      acmat = fullfile(p, [f '.mat']);
      if exist(acmat, 'file')
	s = struct2cell(load(acmat));
	cmap = s{1};
      end
    end
  end
end
if size(cmap, 2)~=3
  warning('Colormap was not an N by 3 matrix')
  cmap = [];
end
