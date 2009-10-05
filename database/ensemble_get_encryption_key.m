function key = ensemble_get_encryption_key(database)
% Returns the encryption key string used by ensemble.
%
% Make sure that the path to the file containing
% the encryption string is defined in ensemble_globals.m. Also
% make sure that you have read permissions to the file containing
% the encryption key.
%
% **********************************************************************************************************  
%
% "Ensemble" is the proprietary property of The Regents of the University of California ("The Regents.")
%
%  Copyright (c) 2005-09 The Regents of the University of California, Davis campus. All Rights Reserved.
%
%  Redistribution and use in source and binary forms, with or without modification, are permitted by
%  nonprofit, research institutions for research use only, provided the conditions in the included
%  license agreement are met.
%
%  Refer to  for the license agreement,
%  
% **********************************************************************************************************
%
% Author(s):
% Sept 30, 2009 - Stefan Tomic, First Version  
  
%ensemble_globals should be a script where
%paths.data_encryption_key is defined
ensemble_globals;

fid = fopen(paths.data_encryption_key,'r');
key = fgetl(fid);
fclose(fid);