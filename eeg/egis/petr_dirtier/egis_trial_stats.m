function [stats, ep] = egis_trial_stats(ses_fname, ep)
%
% [stats, edit_params] = egis_trial_stats(ses_fname, edit_params);

byte_sex = 'ieee-be';

%
%  Read in  information from the session file
%
ses_hdr_offsets_v;

ses_fid = fopen(ses_fname, 'rb', byte_sex);
if ses_fid == -1, error(['Error opening session file' ses_fname]), end

[ses_fhdr,ses_chdr,ses_ename,ses_czeros,ses_cgains,ses_cnames,ses_fcom,ses_ftext, ses_coff]=rd_egis_hdr_v(ses_fid);

%
%  Check for necessary fields in edit_params (ep) structure
%

if nargin < 2
  ep = edit_params;
end

if ~isfield(ep,'start_samp')
  ep.start_samp = 1;
end

if ~isfield(ep,'stop_samp')
  ep.stop_samp = ses_chdr(1,NPoints);
end

if isempty(ep.start_samp), ep.start_samp = 1; end
if isempty(ep.stop_samp), ep.stop_samp = ses_chdr(c, NPoints); end

samp_vect = ep.start_samp:ep.stop_samp;

if ~isfield(ep,'do_eyemove')
  ep.do_eyemove = 0;
end

%
%  Run through all the trial data
%

i = 1;
for c = 1:ses_fhdr(NCells)
  disp(['Processing cell ' int2str(c) ' ...'])
  for t = 1:ses_chdr(c,NTrials)
    disp(['  Trial ' int2str(t) ' ...'])

    stats.cell(i) = c;
    stats.trial(i) = t;
    
    % READ IN DATA FOR SINGLE TRIAL
    trialdata = rd_onetr_allch(ses_fid, ses_coff(c), t, ses_fhdr(NChan), ses_chdr(c, NPoints));

    %
    % CUT OUT THE RELEVANT PORTION
    %
    trialdata = trialdata(samp_vect,:);
    
    % APPLY CALIBRATIONS AND GAINS
    
    if ~any(ses_cgains)
      disp(sprintf('Setting gains to %d', 195))
      ses_cgains = ones(size(ses_cgains)) * 195;
    end
    if ~any(ses_czeros)
      %disp('Setting zeros to zero')
      ses_czeros = zeros(size(ses_czeros));
    end

    stats.cgains = ses_cgains;
    stats.czeros = ses_czeros;

    trialdata = cal_gain(trialdata,ses_cgains,ses_czeros);

    % Compute the mean and standard deviation of the trialdata
    
    stats.vol_mean(i,:) = mean(trialdata);
    stats.vol_std(i,:) = std(trialdata);
    stats.vol_absmax(i,:) = max(max(trialdata, abs(trialdata)));
    stats.vol_diff(i,:) = max(abs(diff(trialdata)));

    % If eye-movement data is desired ...
    
    if ep.do_eyemove
      stats.eye_vert_lft(i) = max(abs(diff(trialdata(:,[ep.suplfteye ...
	    ep.inflfteye]),1,2)));
      stats.eye_vert_rt(i) = max(abs(diff(trialdata(:,[ep.suprteye ...
	    ep.infrteye]),1,2)));
    end % if ep.do_eyemove
    
    i = i + 1;
  end
end

fclose(ses_fid);
