function [cmat, X] = post_process_corrmat(spm_fname,opt)
% plots & returns the design matrix intercorrelation for an SPM analysis
% 
%   [cmat X] = post_process_corrmat(spm_fname)
% 
% REQUIRES
%   opt.figh - handle of figure to print the correlation matrix to
%   opt.names - regressor names (from SPM.xX.name) that you want to
%       correlate. If not provided, all regressors in the design matrix
%       will be correlated
%   opt.plot_corr_txt - text to be printed on figure
%   opt.printfig - (1) print the figure, return values, (0) return values
%   opt.title
%   opt.figstub
%   opt.append_str
% 
% 05/19/06 Petr Janata - modified from a 2004 version of the script.  Added
% support for SPM5 structure of the SPM.mat file.
% 2010.02.16 FB - added opt.names, header documentation

if nargin < 2
  opt = struct();
end

load new_seismic  % colormap

data = load(spm_fname);
if isfield(data,'SPM')
  xX = data.SPM.xX;
  Xnames = xX.name;
  Sess = data.SPM.Sess;
else
  xX = data.xX;
  Xnames = xX.Xnames;
  Sess = data.Sess;
end

% Initialize a figure
if isfield(opt,'figh')
  figure(opt.figh), clf
else
  figure
  opt.figh = gcf;
end

nruns = length(Sess);
for irun = 1:nruns
  % Get the columns corresponding to this run
  cc = Sess(irun).col;
  cr = Sess(irun).row;

  % only plot a subset of correlations?
  if isfield(opt,'names')
    cc = find(ismember(Xnames,opt.names));
  end

  % Get the data corresponding to this run
  X = xX.X(cr,cc);

  % Get the names corresponding to this run
  cnames = Xnames(cc);
  
  % If there is a constant column, remove it
  rem_col = find(sum(X) == size(X,1));
  if ~isempty(rem_col)
    X(:,rem_col) = [];
    cnames{rem_col} = [];
  end
  
  cmat = corrcoef(X);  % calculate the correlation matrix

  % Initialize a subplot
  subplot(nruns,1,irun)
  
  % Plot the matrix
  imagesc(cmat,[-1 1]), axis square

  colormap(new_seismic)
  colorbar

  set(gca,'ytick',1:size(cmat,2))
  set(gca,'yticklabel',cnames)

  if isfield(opt,'plot_corr_txt') && opt.plot_corr_txt
    for ir = 1:size(cmat,1)
      for ic = 1:size(cmat,2)
	text(ic,ir,sprintf('%1.2f', cmat(ir,ic)),'horizontalalign','center')
      end
    end
  end

  set(gca,'xtick',[])
  set(gca,'ticklen',[0 0])

  if nruns == 1
    for ic = 1:size(cmat,2)
      text(ic,size(cmat,1)+0.5,strrep(cnames{ic},'_','\_'),'rotation',-90,'horizontalalign','left')
    end
  end

  if isfield(opt,'title')
    if irun == 1
      title_str = sprintf('%s, Run %d', opt.title, irun);
    else
      title_str = sprintf('Run %d', irun);
    end
    title(title_str,'fontsize',14)
  end
end % for irun

  
try 
  if opt.printfig
    print(sprintf('%s.ps',opt.figstub),'-dpsc','-adobecset', opt.append_str)
    
    %    unix_str = sprintf('ps2pdf %s.ps %s.pdf', opt.figstub, opt.figstub);
    %    unix(unix_str);
    %    unix(sprintf('rm %s.ps >> /dev/null', opt.figstub)); 
  end
catch
end


if nargout > 1
	X = xX.X;
end
