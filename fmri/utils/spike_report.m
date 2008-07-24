function spike_report(info, outdata_root)
% spike_report(info, outdata_root);
%
% Writes spike detection summary data to screen and plots timeseries with
% labelled spikes to a postscript file
%

% 04/07/01 Petr Janata

if nargin < 1
  error('No spike info structure present: run SPIKE_CHECK first')
end

if nargin < 2
  outdata_root = './';
end

fname = fullfile(outdata_root,sprintf('spikereport_%s_%s.ps', datestr(datenum(now),1), datestr(datenum(now),15)));

for isub = 1:length(info)
  if any(info(isub).numbad)
    disp(sprintf('\nSPIKE WARNING: Subject %s', info(isub).id))
    nrun = length(info(isub).run);
    figure(1), clf
    for irun = 1:nrun
      if info(isub).numbad(irun)
	disp(sprintf('\tRun %d: %d spikes', irun, info(isub).numbad(irun)))
      end
      % Plot the data
      subplot(nrun,1,irun)
      plot(info(isub).data{irun})
      bi = info(isub).badidx{irun};
      hold on
      scatter(bi, ones(size(bi))*info(isub).mean(irun)+3*info(isub).std(irun),'*')
      
      if irun == nrun
	ylabel('Signal intensity')
	xlabel('Volume #')
      end
      
      title(sprintf('Subject: %s, Exam: %s, Series: %s, Overall EPI run: %d', info(isub).id, info(isub).exam{irun}, info(isub).series{irun}, irun))
    end
    
    print(gcf,'-dps','-append',fname)
    
  end
end
