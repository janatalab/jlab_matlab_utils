function [onsets,durs] = fmri_regress_timespan(pinfo,minfo,onsets,durs)

% limits onsets and durations by timespan ratings
% 
%   [names,vals] = fmri_regress_timespan(pinfo,minfo,ons,durs)
% 
% takes presentation data, model information, onsets and durations, finds
% and calculates timespan ratings, and then adjusts the onsets and durs to
% only reflect the timespan that was rated for the given response. See
% also: fmri_regress_stim.m
% 
% FB 2010.01.27

  if isfield(minfo,'timespan_defs')
    tsd = minfo.timespan_defs;
  else
    error('no timspan defs found');
  end
  
  pc = set_var_col_const(pinfo.vars);
  
  % calculate timespan ratings from presentation data
  tsfind = strfind(pinfo.data{pc.EVENT_CODE},'timespan_cue');
  tsidxs = find(~cellfun(@isempty,tsfind));
  
  for j=1:length(tsidxs);
    codes = pinfo.data{pc.RESP_CODE}{tsidxs(j)};
    pos = tsd.start; lvl = 0; start = 1; stop = tsd.nbins;
    for k=1:length(codes)
      code = str2num(codes{k});
      switch code
        case tsd.resp.enter
          if lvl == 0
            if pos == tsd.nbins
              pos = tsd.nbins - 1;
            end
            start = pos;
            pos = pos + 1;
            lvl = lvl+1;
          elseif lvl == 1
            stop = pos;
            lvl = lvl+1;
          end % if lvl == 0
        case {tsd.resp.left,tsd.resp.right}
          if lvl == 2
            if code == tsd.resp.left
              pos = tsd.start;
              start = 1;
              stop = tsd.nbins;
              lvl = 0;
            elseif code == tsd.resp.right
              lvl = lvl + 1;
              break
            end
          elseif lvl < 2
            if code == tsd.resp.left
              pos = max(pos - 1,start);
            elseif code == tsd.resp.right
              pos = min(pos + 1,tsd.nbins);
            end
          end % if lvl == 2
      end % switch code
    end % for k=1:length(codes
    start = start/tsd.nbins;
    stop = stop/tsd.nbins;
    
    % link this timespan judgment with a stimulus
    time = pinfo.data{pc.RUN_REL_TIME}(tsidxs(j))/1000;
    oidx = find(onsets < time,1,'last');
    
    % modify the onset/dur for this stimulus
    onsets(oidx) = onsets(oidx)+minfo.music_dur*start;
    durs(oidx) = minfo.music_dur*(stop-start);
  end % for j=1:length(tsidxs
  