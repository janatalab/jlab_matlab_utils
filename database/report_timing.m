function rp = report_timing(as,rp)
% Generates timing reports for given datasets.
% report_timing.m
%
% Generate various reports and plots pertaining to timing aspects of the data
% contained in the analysis structure (as).  See repack_formdata.m for more
% information on the as structure.
%
% rp is a structure containing the report parameters, including the id of the
% report type. rp can be a structure array in which each element corresponds to
% a different report type.
%
% Current report types:
%
%   'time_to_completion'
%

% 09/12/06 Petr Janata - modified from report_autobio

nrep = length(rp);

for irep = 1:nrep
  rep_id = rp(irep).id;
  switch rep_id
    case 'time_to_completion'
      starttime = [];
      stoptime = [];
      
      % Generate a matrix with subjects in columns
      if isfield(as,'nnum')
	if isfield(as.num.datenum,'by_question')
	  data_dims = size(as.num.datenum.by_question);
	  data = reshape(as.num.datenum.by_question, [data_dims(1) ...
		prod(data_dims(2:3))])';
	  
	  exist_mask = data > 0;
	  data(~exist_mask) = NaN;
	  
	  starttime = nanmin(data);
	  stoptime = nanmax(data);
	end
      elseif isfield(as,'sess')
	starttime = as.sess.start_time;
	stoptime = as.sess.stop_time;
      end
      
      % Generate a figure showing the distribution of completion times
      elapsed_time = etime(datevec(stoptime),datevec(starttime));
      min_time = min(elapsed_time);
      max_time = max(elapsed_time);
      mean_time = mean(elapsed_time);
      median_time = median(elapsed_time);
      
      fprintf('Distribution of times in raw data:\n');
      fprintf('Minimum elapsed time (min): %1.2f\n', min_time/60);
      fprintf('Maximum elapsed time (min): %1.2f\n', max_time/60);
      fprintf('Mean elapsed time(min): %1.2f\n', mean_time/60);
      fprintf('Median elapsed time(min): %1.2f\n', median_time/60);
      
      % Limit the display by range
      if isfield(rp,'range')
	hist_min = rp.range(1);
	hist_max = rp.range(2);
      else
	hist_min = min_time;
	hist_max = max_time;
      end
      
      rp.results.starttime = starttime;
      rp.results.stoptime = stoptime;
      rp.results.elapsed_time = elapsed_time;
      
      % Generate a mask of for data are in range
      valid_mask = (elapsed_time >= hist_min) & (elapsed_time <= hist_max);
      
      xscale = linspace(hist_min,hist_max,100);
      [histdata] = histc(elapsed_time(valid_mask), xscale);
      
      if rp.plot_data
	figure(1),clf
	bar(xscale/60,histdata);
	if isfield(rp,'tick_interval')
	  set(gca,'xtick',min(xscale):rp.tick_interval:max(xscale));
	end
      end
      
    otherwise
      fprintf('No handling for report type: %s\n', rep_id);
  end
end % for irep=
