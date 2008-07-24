function data = batch_shuffle_check(dataroot, sinfo, outfstub, epidir)
% batch_shuffle_check(dataroot, sinfo)
%
% Performs slice-shuffle checking for all subjects specified in sinfo.
%
% For each functional run, the script goes through all volumes slice by slice
% and creates a color image in which the columns represent time and the rows
% represent slices. The color of each slice at each time point relects the
% summed image intensity for that slice.

% Input:
%   dataroot -- directory in which all subject directories live

%   sinfo -- a structure array in which each element contains information about
%            the subject. Required fields are:
%     .id -- the subject string, e.g. 'sub01', which is the name of the
%            directory in which the subject's data live.
%     .use_runs -- a list of EPI runs to analyze
%  
%   outfstub -- the root name of the file into which the images will be
%               written.  Assuming that ps2pdf works, a PDF file is created
%               that contains the data for all subjects.
%
%   epidir -- the name of the subdirectory within subject directories that the
%             data exist in.  Default is 'epi'.  It is assumed that individual
%             run directories are called run? where ? is the run number.
%
% Output:
%   data -- the actual data

% 09/27/03 Petr Janata

if nargin < 4
  epidir = 'epi';
end

if nargin < 3
  outfstub = 'shuffle_check';
end

PRINT_TO_FILE = 1;

proc_sub = 1:length(sinfo);
nsub = length(proc_sub);

for isub = 1:nsub
  sub_idx = proc_sub(isub);
  
  runs = sinfo(sub_idx).use_runs;
  nruns = length(runs);

  figure(1), clf
  
  for irun = 1:nruns
    run_idx = runs(irun);
    epi_dir = fullfile(dataroot,sinfo(sub_idx).id,epidir,sprintf('run%d', run_idx));
    P = spm_get('Files',epi_dir,'*.img');
    nvol = length(P);
  
    % Get the number of slices from the first volume
    V = spm_vol(deblank(P(1,:)));
    nslice = V.dim(3);
    
    data{sub_idx,irun} = zeros(nslice,nvol);
    for ivol = 1:nvol
      fprintf('Subject: %s, Run %d, Vol %d\n', sinfo(sub_idx).id,run_idx, ivol)
      V = spm_vol(deblank(P(ivol,:)));
      for islice = 1:nslice
	tmp = spm_slice_vol(V,spm_matrix([0 0 islice]),V.dim(1:2),0);
	data{sub_idx,irun}(islice,ivol) = sum(tmp(:));
      end % for islice=
    end % for ivol=
    
    subplot(nruns,1,irun);
    imagesc(data{sub_idx,irun})
    set(gca,'xtick',[0:10:nvol],'xticklabel','')
    
    tick_locs = get(gca,'xtick');
    use_ticks = 1:5:length(tick_locs);
    for itick = 1:length(use_ticks)
      tick_val = tick_locs(use_ticks(itick));
      text(tick_val,nslice+1,num2str(tick_val),'rotation',-90)
    end
    
    grid
    title(sprintf('Subject: %s, Run %d', sinfo(sub_idx).id,run_idx))
    
  end % for irun=

  if PRINT_TO_FILE
    if isub ~= 1
      append_str = '-append';
    else
      append_str = '';
    end
    print('-dpsc',append_str,sprintf('%s.ps',outfstub))
  end
end % for isub =

if PRINT_TO_FILE
  unix_str = sprintf('ps2pdf %s.ps %s.pdf', outfstub,outfstub);
  unix(unix_str);
end
