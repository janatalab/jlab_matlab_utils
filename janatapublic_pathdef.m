function p = janatapublic_pathdef
%PATHDEF Search path defaults.
%   PATHDEF returns a string that can be used as input to MATLABPATH
%   in order to set the path.

  
%   Copyright 1984-2000 The MathWorks, Inc.
%   $Revision: 1.4 $ $Date: 2000/06/01 16:19:21 $

% PATH DEFINED HERE -- Don't change this line.

switch computer
  case {'PCWIN','PCWIN64'}
    path_sep = ';';
  otherwise
    path_sep = ':';
end

this_path = which('janatapublic_pathdef.m');
public_path = this_path(1:findstr(this_path,'janatapublic_pathdef.m')-1);

p = [...
      fullfile(public_path,['utils' path_sep]),...
      fullfile(public_path,'utils',['parser' path_sep]),...
      fullfile(public_path,'utils',['physio' path_sep]),...
      fullfile(public_path,'utils',['signal' path_sep]),...
      fullfile(public_path,'utils',['plot' path_sep]),...
      fullfile(public_path,['database' path_sep]), ...
... %      fullfile(public_path,'database', ['mysql_cpp' path_sep]), ...
      fullfile(public_path,'database',['physio' path_sep]),...
      fullfile(public_path,'fmri',['fsl' path_sep]),...
      fullfile(public_path,'fmri',['spm5_struct_templates' path_sep]),...
      fullfile(public_path,'fmri',['utils' path_sep]),...
      fullfile(public_path,'fmri','utils',['physio' path_sep]),...
      fullfile(public_path,'fmri','utils',['GE2SPM' path_sep]),...
      fullfile(public_path,'fmri','utils',['vanhorn' path_sep]),...
      fullfile(public_path,['colormaps' path_sep]),...
      fullfile(public_path,'eeg',['eeglab_petrmods' path_sep]),...
      fullfile(public_path,'eeg',['fixes' path_sep]),...
      fullfile(public_path,'eeg',['utils' path_sep]),...
      ];
