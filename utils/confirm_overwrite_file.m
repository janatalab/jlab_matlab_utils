function confirmation = confirm_overwrite_file(filepath)
%
% confirmation = confirm_overwrite_file(filepath)
%
% confirmation is 1 or 0, indicating whether or not filepath exists after
% this function is run. (in other words, if the file was deleted by
% this function or the file 
%
% checks to see if a file designated by filepath exists. If it does
% exist, prompts the user whether he/she want's to overwrite the
% file. If the user answers 'y' the file is deleted. This utility
% is useful when a function makes repeated 'print' or fprintf
% calls, appending to a plot or log file and it is desireable to
% delete the plot or log file at the begninning of the
% process. User confirmation ensures that useful plots or logs are
% not simply overwritten,
%
% 12 July 2007, Stefan Tomic

if(exist(filepath,'file'))

 confirmOverwrite = '';
  while(~ismember(confirmOverwrite,{'y','n'}))
    confirmOverwrite = input(sprintf('File %s exists, overwrite? (y/n): ',filepath),'s');
  end
  
  if(strcmp(confirmOverwrite,'y'))
    delete(filepath);
    confirmation = 1;
  else
    confirmation = 0;
  end

else 
  confirmation = 1;
end






return
