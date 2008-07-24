function bad_stim_ids = check_ensemble_mp3s
% Checks mp3s to make sure they are valid mp3 files.
% needs cleaning!!
  bad_stim_ids = [];
  stim_path = '/var/www/html/questionnaire/stimuli';
  
  mysql_make_conn;
  sql_get_stims = sprintf('select stimulus_id,artist,name,location from stimulus where file_format = ''mp3''');
  [stimulus_id,artist,name,location] = mysql(sql_get_stims);
  
  nstims = length(stimulus_id);
  
  for istim = 1:nstims
    stim_fullname = fullfile(stim_path,char(location(istim))); 
    eval(sprintf('!mp3info -x %s;',stim_fullname),'bad_stim_ids = [bad_stim_ids;stimulus_id(istim)];');
   
  end
  
  
  mysql('close');
  
  
end
