% rfchop_start.m

global dataroot proc_subs

SCAN_OFFSET = 0;

CONVERT_FORMAT = []; % 
GLOBAL_SPECT = 1;
DO_CALC = 1;
PLOT_SPECT = 1;
LOWPASS = 1;

mask_threshold = 500;
Fs = 1/2;

proc_subs = 1:2;

s = get_rfchop_sinfo;

%
% Convert data from GE format into AVW
%

if CONVERT_FORMAT
  convert_format(s(CONVERT_FORMAT),dataroot,dataroot, SCAN_OFFSET);
end

if GLOBAL_SPECT
  if DO_CALC
    clear avgmag stdmag ss
  end
  for isub = 1:length(proc_subs)
    subidx = proc_subs(isub);

    fprintf('Subject: %s\n', s(subidx).id);
    
    subdir = fullfile(dataroot, s(subidx).id);
    nruns = length(s(subidx).use_runs);
    
    if DO_CALC
      for irun = 1:nruns
	runidx = s(subidx).use_runs(irun);
	rundir = fullfile(subdir,'epi',sprintf('run%d',runidx));
	
	fprintf('Run (%d/%d): %d\n', irun, nruns, runidx)
	
	% Load all of the images
	flist = dir(fullfile(rundir,sprintf('%s*i0*.img',s(subidx).id)));
	V = spm_vol(fullfile(rundir,flist(1).name));
	if length(flist) ~= s(subidx).nvol(runidx)
	  error('Wrong number of files in directory')
	end
	
	fprintf('Loading data\n')
	data = zeros(s(subidx).nvol(runidx),prod(V.dim(1:3)));
	nvol = s(subidx).nvol(runidx);
	for ivol = 1:nvol
	  fprintf('\t%s\n', flist(ivol).name);
	  fid = fopen(fullfile(rundir,flist(ivol).name),'rb');
	  data(ivol,:) = fread(fid,spm_type(V.dim(4)))';
	  fclose(fid);
	end % for ivol=
	
	% Get voxels that exceed threshold
	vox = find(data(1,:) > mask_threshold);

	% Remove mean from each voxel
	data = detrend(data,'constant');
	
	% Filter if desired
	if LOWPASS
	  [b,a] = butter(5,0.8);
	  for ivox = 1:length(vox)
	    data(:,vox(ivox)) = filtfilt(b,a,data(:,vox(ivox)));
	  end
	end
	
	% Compute FFTs
	fprintf('Computing spectra ...\n');
	magdata = abs(fft(data(:,vox)))/nvol;
	avgmag{isub}(:,irun) = mean(magdata,2);
	stdmag{isub}(:,irun) = std(magdata,[],2);

	
	% Get SS across all voxels
	ss{isub}(irun,:) = sum(data.^2);
	avg_ss{isub}(irun) = mean(ss{isub}(irun,:),2);
	std_ss{isub}(irun) = std(ss{isub}(irun,:),[],2);
      end % for irun=
    end % if DO_CALC
    
    if PLOT_SPECT
      figure(isub),clf

      nvol = s(subidx).nvol(runidx);

      fscale = (0:fix(nvol/2)-1)/fix(nvol/2)*(Fs/2);
      plot(fscale,(avgmag{isub}(1:fix(nvol/2),:).^2)./repmat(avg_ss{isub},fix(nvol/2),1))
      xlabel('Frequency (Hz)')
      ylabel('Power')
      legend(strrep(cond_names(s(subidx).cond_order(s(subidx).use_runs)),'_','\_'))
    end
  end % for isub=
end % if GLOBAL_SPECT

