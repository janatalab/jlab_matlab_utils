% timecourse.m
%
% Extracts time-course from an ROI
%

global BCH; %- used as a flag to know if we are in batch mode or not.

USE_INTERFACE = 1;
USE_YMAD = 0;

rootpath = '/home03/handy/imagery/';
spmpath = 'analysis';
sub_ids = ...
    {'02jun00CW', ...
     '03aug00GG', ...
     '10aug00EH', ...
     '11aug00WS', ...
     '14jun00OMG', ...
     '15jun00JE', ...
     '15jun00SN', ...
     '17jul00RC', ...
     '18jul00CS'}; 

series_map.episodic = [ ...
    3 6; ...
    2 5; ...
    2 6; ...
    3 6; ...
    3 6; ...
    2 5; ...
    3 6; ...
    2 5; ...
    3 6];

series_map.semantic = [ ...
    1 4; ...
    3 6; ...
    3 6; ...
    1 4; ...
    1 4; ...
    3 6; ...
    1 4; ...
    3 6; ...
    1 4];

task_vols = [16:34 50:68 84:102 118:136];
rest_vols = [1:15 35:49 69:83 103:117];

nser = 6;
nsub = length(sub_ids);


if USE_INTERFACE

  %
  %  Deal with getting the cluster of interest from the image of interest.  This
  %  is the cluster that will be used to pull data out of individual subject files.
  %

  % from spm_regions.m

  Finter = spm_figure('GetWin','Interactive');
  Fgraph = spm_figure('GetWin','Graphics');
  set(Finter,'Name','VOI time-series extraction')

  %-Find nearest voxel [Euclidean distance] in point list with data saved
  % in Y.mad, and update GUI location
  %-----------------------------------------------------------------------
  if ~length(SPM.XYZmm)
    spm('alert!','No suprathreshold voxels!',mfilename,0);
    Y = []; xY = [];
    return
  elseif exist(fullfile(SPM.swd,'Y.mad')) ~= 2
    spm('alert!','No raw data saved with this analysis!',mfilename,0);
    Y = []; xY = [];
    return
  elseif ~any(SPM.QQ)
    spm('alert!','No raw data saved for any suprathreshold location!',...
	mfilename,0);
    Y = []; xY = [];
    return
  end

  % from spm_VOI.m and spm_regions.m
  xyzmm   = spm_results_ui('GetCoords');

  [xyzmm,i] = spm_XYZreg('NearestXYZ',xyzmm,SPM.XYZmm);
  spm_results_ui('SetCoords',xyzmm);
  A     = spm_clusters(SPM.XYZ);
  Q     = find(A == A(i));
  str   = sprintf('%0.0f voxel cluster',length(Q))
  D     = SPM.XYZ(:,Q);
  S     = length(Q);
  
  q       = find(SPM.QQ(Q));

  if any(SPM.QQ(Q)==0)
    spm('alert"',{...
        sprintf('Don''t have raw data for all %d suprathreshold',length(Q)),...
	sprintf('voxels.'), ...
	sprintf('Proceeding using the %d voxels that are.',length(q)),...
		 },mfilename,sqrt(-1));
  end
else
  D = global_cluster;
end % if USE_INTERFACE

%
% Read the data from one of two sources.  
%   One is the Y.mad file
%   The other is from raw data files

%-Get (approximate) raw data y from Y.mad file
%-NB: Data in Y.mad file is compressed, and therefore not fully accurate
%-----------------------------------------------------------------------

clear y

if USE_YMAD
  y       = spm_extract(fullfile(SPM.swd,'Y.mad'),SPM.QQ(Q(q)));
else

  for isub = 1:nsub
    clear tmp
    for iser = 1:nser
      switch sub_ids{isub}
       case '14jun00OMG'
	tmp{iser} = spm_get('Files', ... 
			    sprintf('%s%s/spm/functional/series%02d/', ...
				    rootpath, sub_ids{isub}, iser), ...
			    'sns*14jun000MG*.img');
       otherwise
	tmp{iser} = spm_get('Files', ... 
			    sprintf('%s%s/spm/functional/series%02d/', rootpath, sub_ids{isub}, iser), ... 
			    sprintf('sns*%s*.img', sub_ids{isub}));
      end
    end
    filelist{isub} = char(tmp);

%    clear VY
%    fname = fullfile(rootpath, sub_ids{isub}, spmpath, 'SPM.mat');
%    load(fname)
    
%    filelist{isub} = strvcat(VY.fname);
    nvol = size(filelist{isub},1);    
    disp(sprintf('Loading data for subject: %s', sub_ids{isub}))
    for ivol = 1:nvol
      %
      %  Loop through all of the raw data files, grabbing the relevant voxels
      %
      vol = (spm_vol(deblank(filelist{isub}(ivol,:))));
      y(:,ivol,isub) = spm_sample_vol(vol, D(1,:),D(2,:),D(3,:),0);
    end
  end  % for isub=
end