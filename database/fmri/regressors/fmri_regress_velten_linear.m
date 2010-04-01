function [names,vals] = fmri_regress_velten_linear(pinfo,minfo,sess)

% generates linear blocked fmri data regressors for the velten mood induction task
% 
%   [names,vals] = fmri_regress_velten_linear(pinfo,minfo,sess)
% 
% called by fmri_generate_regress
% 
% REQUIRES
% 
% RETURNS
%   names = cell array of six regressor names
%   vals = volume X regressor matrix containing motion regressors
% 
% FB 2010.03.16

regid = pinfo.regid;

% init output vars
names = {};
vals = [];

rnum = pinfo.irun;
if isempty(strfind(pinfo.presfname,'velten')), return, end
vidx = strfind(pinfo.presfname,'velten')+7;
eidx = strfind(pinfo.presfname(vidx:end),'_');
emo  = pinfo.presfname(vidx:vidx+eidx-2);

pc = set_var_col_const(pinfo.vars);

% get all event blocks
rfilt.include.all.EVENT_CODE = {'velten_.*_block.*'};
rinfo = ensemble_filter(pinfo,rfilt);
nr = length(rinfo.data{1});

% get xhair blocks
xfilt.include.all.EVENT_CODE = {'xhair'};
xinfo = ensemble_filter(pinfo,xfilt);
nx = length(xinfo.data{1});

if (nr == 0)
    % check to see if stim blocks are not labeled
    blfilt.include.all.EVENT_CODE = {'block_*'};
    blfilt.include.all.EVENT_TYPE = {'Picture'};
    bldata = ensemble_filter(pinfo,blfilt);
    nbl = length(bldata.data{1});

    picfilt.include.all.EVENT_CODE = {'pic_*'};
    picfilt.include.all.EVENT_TYPE = {'Picture'};
    picdata = ensemble_filter(pinfo,picfilt);
    npic = length(picdata.data{1});

    if ~nbl || ~npic
      fprintf(1,'mood trend regressor skipped: no onsets found\n');
      return
    elseif npic < nbl
      fprintf(1,['fewer stim (%d) than block onsets (%d) for '...
          'regressor %s\n'],npic,nbl,regid);
      return
    elseif nbl <= nx
      fprintf(1,['fewer or equal blocks (%d) and xhairs (%d) '...
          'for regressor (%s)\n'],nbl,nx,regid);
      return
    end

    % there are blocks, and there are pictures, line them up
    % find the first pic after every block, this is your stim block
    % onset. assume the last blon is a final rest block, and is not
    % followed by a stim block.
    blons = bldata.data{pc.RUN_REL_TIME}/1000;
    picons = picdata.data{pc.RUN_REL_TIME}/1000;
    xhons = xinfo.data{pc.RUN_REL_TIME}/1000;

    amp = ones(1,nbl-1);
    onsets = [];
    offsets = [];

    for ib=1:nbl-1
      lons = find(picons > blons(ib));
      if picons(lons(1)) < xhons(ib)
        onsets = [onsets picons(lons(1))];
        offsets = [offsets xhons(ib)];
      else
        fprintf(1,['onset (%ds) later than offset (%ds) for '...
            'block %d, regressor %s\n'],picons(lons(1)),...
            xhons(ib),ib,regid);
      end
    end
else
    % number of blocks
    if (nx == 0)
        fprintf('no crosshairs found, skipping %s\n',regid);
        return
    elseif (nr == nx)
        nb = nr;
    elseif ((nx > nr) && ~isempty(min(rinfo.data{pc.RUN_REL_TIME}) > ...
        min(xinfo.data{pc.RUN_REL_TIME})))
        % first xhair is before first block onset ... delete all
        % xhairs before the first block onset
        while (min(rinfo.data{pc.RUN_REL_TIME}) > ...
          min(xinfo.data{pc.RUN_REL_TIME}))
            for iv = 1:length(xinfo.data)
              xinfo.data{iv}(1) = [];
            end
        end
        nb = nr;
    else
        nb = min(nr,nx);
    end

    % build regressor, if nb > 0 and if rnum <= nb
    if (nb == 0)
        return
    end
            
    % get onsets and offsets
    onsets = rinfo.data{pc.RUN_REL_TIME}/1000;
    offsets = xinfo.data{pc.RUN_REL_TIME}/1000;

    % if more offsets or onsets than blocks, throw a warning, but
    % assume for the moment that there are extra items at the end,
    % and cut them out
    if nx > nb
      offsets = offsets(1:nb);
      fprintf(1,'more offsets (%d) than blocks (%d) in %s\n',...
          nx,nb,regid);
    end
    if nr > nb
      onsets = onsets(1:nb);
      fprintf(1,'more onsets (%d) than blocks (%d) in %s\n',...
          nr,nb,regid);
    end
end

% find out what emotions, and how many blocks there are
nons = length(onsets);
names = cell(1,nons);
vals = zeros(pinfo.scanner.actual_nvol,nons);
TR = pinfo.scanner.TR;

for j=1:nons
  names{j} = sprintf('velten_%s_lin%d',emo,j);
  ontr  = round(onsets(j)/TR);
  offtr = round(offsets(j)/TR)-1;
  nTRs  = round(offtr - ontr) + 1;
  steps = 1/(nTRs - 1);
  vals(ontr:offtr,j) = [0:steps:1];
end
