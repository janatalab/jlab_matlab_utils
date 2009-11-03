function [cmat, X] = post_process_corrmat_fsl(fsldir,opt)
% plots design matrix and regressor covariance structure for FSL analyses
% 
%   [cmat X] = post_process_corrmat_fsl(fsldir)
% 
% REQUIRES
%   fsldir - path to fsl analysis directory, within which should exist:
%       design.mat, design.fsf, fsf.mat
%   opt.printfig - flag, 1 or 0, to print or not print the figure to file
%   opt.figstub - file name stub (w/o extension) for figure filenames
%   opt.appendstr - cell array of strings of additional settings to send to
%       print when printing the figures to file
%   opt.title - figure title stub (to which will be appended 'design
%       matrix' or 'regressor correlations'
%   opt.plot_corr_txt - plot text labels on regressor correlations?
% 
% RETURNS
%   cmat - correlation matrix for design matrix regressors
%   X - design matrix
% 
% 05/19/06 Petr Janata - modified from a 2004 version of the script.  Added
% support for SPM5 structure of the SPM.mat file.
% 10/28/09 Fred Barrett - adapted from post_process_corrmat (SPM only) for
% use with FSL

if nargin < 2
  opt = struct();
end

load new_seismic  % colormap

if ~exist(fsldir,'dir')
  error('fsldir (%s) not found',fsldir);
end

matfname = fullfile(fsldir,'fsf.mat');
fsffname = fullfile(fsldir,'design.fsf');
desfname = fullfile(fsldir,'design.mat');

% get FEAT fsf structure
if ~exist(matfname,'file') && ~exist(fsffname,'file')
  error('FSL fsf FEAT structure not found in %s',fsldir);
elseif exist(matfname,'file')
  load(matfname);
elseif exist(fsffname,'file')
  error(['reading straight from design.fsf not yet supported ...\n' ...
      'please provide an fsf.mat file']);
  % otherwise, you may load the fsffname file here, and parse the file to
  % get an fsf structure    
end

evnames = {fsf.ev.name};

% get design matrix
if exist(desfname,'file')
  desmat = loadtxt(desfname,'skipline',5);
  X = cell2mat(desmat);
else
  error('could not find design matrix (%s)',desfname);
end

% % % % % % PLOT: design matrix

% Initialize a figure
figure

% Plot the design matrix
imagesc(X);

% colormap(new_seismic)
% colorbar

set(gca,'xtick',[])
set(gca,'ticklen',[0 0])

for ic = 1:size(X,2)
  text(ic,size(X,1)+0.5,strrep(evnames{ic},'_','\_'),'rotation',-90,'horizontalalign','left')
end

% add title
if isfield(opt,'title')
  title(sprintf('%s design matrix',opt.title),'fontsize',14)
else
  title('design matrix','fontsize',14);
end

% %  print to file

try 
  if opt.printfig
    print(sprintf('%s_desmat.ps',opt.figstub),'-dpsc','-adobecset', opt.append_str)
    
    %    unix_str = sprintf('ps2pdf %s.ps %s.pdf', opt.figstub, opt.figstub);
    %    unix(unix_str);
    %    unix(sprintf('rm %s.ps >> /dev/null', opt.figstub)); 
  end
catch
end

% % % % % % PLOT: regressor correlations

% Initialize a figure
figure

% If there is a constant column, remove it
rem_col = find(sum(X) == size(X,1));
if ~isempty(rem_col)
    X(:,rem_col) = [];
    evnames(rem_col) = [];
end
  
cmat = corrcoef(X);  % calculate the correlation matrix

% Plot the matrix
imagesc(cmat,[-1 1]), axis square

% colormap(new_seismic)
colorbar

set(gca,'ytick',1:size(cmat,2))
set(gca,'yticklabel',evnames)

if isfield(opt,'plot_corr_txt') && opt.plot_corr_txt
  for ir = 1:size(cmat,1)
    for ic = 1:size(cmat,2)
      text(ic,ir,sprintf('%1.2f', cmat(ir,ic)),'horizontalalign','center')
    end
  end
end

set(gca,'xtick',[])
set(gca,'ticklen',[0 0])

for ic = 1:size(cmat,2)
  text(ic,size(cmat,1)+0.5,strrep(evnames{ic},'_','\_'),'rotation',-90,'horizontalalign','left')
end

if isfield(opt,'title')
  title(sprintf('%s regressor correlations',opt.title),'fontsize',14)
else
  title('design matrix','fontsize',14);
end

% %  print to file

try
  if opt.printfig
    print(sprintf('%s_desmat_corr.ps',opt.figstub),'-dpsc','-adobecset', opt.append_str)
    
    %    unix_str = sprintf('ps2pdf %s.ps %s.pdf', opt.figstub, opt.figstub);
    %    unix(unix_str);
    %    unix(sprintf('rm %s.ps >> /dev/null', opt.figstub)); 
  end
catch
end
