% ensemble_fmri_set_sess_paths

% FB 2008.08.29 - adapted from PJ/autobio_fmri version
% adding directories generated to outdata.data{paths_idx}
% FB 2009.03.11 - using ensemble_add_data_struct_row

%
% Deal with directory paths
%

if exist('pathdata','var') && ~isempty(strmatch('sess_indir',...
        pathdata.data{pacol.path_type}))
  sess_indir = pathdata.data{pacol.path}{strmatch('sess_indir',...
      pathdata.data{pacol.path_type})};
else
  sess_indir = fullfile(sub_indir,session_stub);
  check_dir(sess_indir,1);
  
  if exist('paths_idx','var')
    outdata.data{paths_idx} = ensemble_add_data_struct_row(...
        outdata.data{paths_idx},'subject_id',subid,'session',isess,...
        'run',0,'path_type','sess_indir','path',sess_indir);
  end
end

if exist('pathdata','var') && ~isempty(strmatch('sess_outdir',...
        pathdata.data{pacol.path_type}))
  sess_outdir = pathdata.data{pacol.path}{strmatch('sess_outdir',...
      pathdata.data{pacol.path_type})};
else
  sess_outdir = fullfile(sub_outdir,session_stub);
  check_dir(sess_outdir,1);
  
  if exist('paths_idx','var')
    outdata.data{paths_idx} = ensemble_add_data_struct_row(...
        outdata.data{paths_idx},'subject_id',subid,'session',isess,...
        'run',0,'path_type','sess_outdir','path',sess_outdir);
  end
end

if exist('pathdata','var') && ~isempty(strmatch('behav_indir',...
        pathdata.data{pacol.path_type}))
  behav_indir = pathdata.data{pacol.path}{strmatch('behav_indir',...
      pathdata.data{pacol.path_type})};
else
  % For the behavior files
  behav_indir = fullfile(sess_indir,'behavior');
  check_dir(behav_indir);

  if exist('paths_idx','var')
    outdata.data{paths_idx} = ensemble_add_data_struct_row(...
        outdata.data{paths_idx},'subject_id',subid,'session',isess,...
        'run',0,'path_type','behav_indir','path',behav_indir);
  end
end

if exist('pathdata','var') && ~isempty(strmatch('behav_outdir',...
        pathdata.data{pacol.path_type}))
  behav_outdir = pathdata.data{pacol.path}{strmatch('behav_outdir',...
      pathdata.data{pacol.path_type})};
else
  % For the behavior files
  behav_outdir = fullfile(sess_outdir,'behavior');
  check_dir(behav_outdir);

  if exist('paths_idx','var')
    outdata.data{paths_idx} = ensemble_add_data_struct_row(...
        outdata.data{paths_idx},'subject_id',subid,'session',isess,...
        'run',0,'path_type','behav_outdir','path',behav_outdir);
  end
end

if exist('pathdata','var') && ~isempty(strmatch('physio_indir',...
        pathdata.data{pacol.path_type}))
  physio_indir = pathdata.data{pacol.path}{strmatch('physio_indir',...
      pathdata.data{pacol.path_type})};
else
  % For the physio files
  physio_indir = fullfile(sess_indir,'physio');
  check_dir(physio_indir);

  if exist('paths_idx','var')
    outdata.data{paths_idx} = ensemble_add_data_struct_row(...
        outdata.data{paths_idx},'subject_id',subid,'session',isess,...
        'run',0,'path_type','physio_indir','path',physio_indir);
  end
end

if exist('pathdata','var') && ~isempty(strmatch('physio_outdir',...
        pathdata.data{pacol.path_type}))
  physio_outdir = pathdata.data{pacol.path}{strmatch('physio_outdir',...
      pathdata.data{pacol.path_type})};
else
  % For the physio files
  physio_outdir = fullfile(sess_indir,'physio');
  check_dir(physio_outdir);

  if exist('paths_idx','var')
    outdata.data{paths_idx} = ensemble_add_data_struct_row(...
        outdata.data{paths_idx},'subject_id',subid,'session',isess,...
        'run',0,'path_type','physio_outdir','path',physio_outdir);
  end
end

if exist('pathdata','var') && ~isempty(strmatch('epi_indir',...
        pathdata.data{pacol.path_type}))
  epi_indir = pathdata.data{pacol.path}{strmatch('epi_indir',...
      pathdata.data{pacol.path_type})};
else
  % For the EPI files
  % Here I have to do something a bit cludgy, because my original and
  % analyzed data branch at the session level.
  epi_indir = fullfile(sess_indir, 'orig');
  check_dir(epi_indir,1);

  if exist('paths_idx','var')
    outdata.data{paths_idx} = ensemble_add_data_struct_row(...
        outdata.data{paths_idx},'subject_id',subid,'session',isess,...
        'run',0,'path_type','epi_indir','path',epi_indir);
  end
end

if exist('pathdata','var') && strmatch('epi_outdir',...
        pathdata.data{pacol.path_type})
  epi_outdir = pathdata.data{pacol.path}{strmatch('epi_outdir',...
      pathdata.data{pacol.path_type})};
else
  epi_outdir = fullfile(sess_outdir, 'epi');
  check_dir(epi_outdir,1);

  if exist('paths_idx','var')
    outdata.data{paths_idx} = ensemble_add_data_struct_row(...
        outdata.data{paths_idx},'subject_id',subid,'session',isess,...
        'run',0,'path_type','epi_outdir','path',epi_outdir);
  end
end

if exist('pathdata','var') && strmatch('anat_indir',...
        pathdata.data{pacol.path_type})
  anat_indir = pathdata.data{pacol.path}{strmatch('anat_indir',...
      pathdata.data{pacol.path_type})};
else
  anat_indir = fullfile(sess_indir, 'orig');

  if exist('paths_idx','var')
    outdata.data{paths_idx} = ensemble_add_data_struct_row(...
        outdata.data{paths_idx},'subject_id',subid,'session',isess,...
        'run',0,'path_type','anat_indir','path',anat_indir);
  end
end

if exist('pathdata','var') && strmatch('anat_outdir',...
        pathdata.data{pacol.path_type})
  anat_outdir = pathdata.data{pacol.path}{strmatch('anat_outdir',...
      pathdata.data{pacol.path_type})};
else
  anat_outdir = fullfile(sess_outdir, 'anatomy');
  check_dir(anat_outdir,1);

  if exist('paths_idx','var')
    outdata.data{paths_idx} = ensemble_add_data_struct_row(...
        outdata.data{paths_idx},'subject_id',subid,'session',isess,...
        'run',0,'path_type','anat_outdir','path',anat_outdir);
  end
end

if exist('pathdata','var') && strmatch('anal_outdir',...
        pathdata.data{pacol.path_type})
  anal_outdir = pathdata.data{pacol.path}{strmatch('anal_outdir',...
      pathdata.data{pacol.path_type})};
else
  anal_outdir = fullfile(sess_outdir, 'analyses');
  check_dir(anal_outdir);

  if exist('paths_idx','var')
    outdata.data{paths_idx} = ensemble_add_data_struct_row(...
        outdata.data{paths_idx},'subject_id',subid,'session',isess,...
        'run',0,'path_type','anal_outdir','path',anal_outdir);
  end
end

if exist('pathdata','var') && strmatch('spm_outdir',...
        pathdata.data{pacol.path_type})
  spm_outdir = pathdata.data{pacol.path}{strmatch('spm_outdir',...
      pathdata.data{pacol.path_type})};
else
  spm_outdir = fullfile(anal_outdir, 'spm');
  check_dir(spm_outdir);

  if exist('paths_idx','var')
    outdata.data{paths_idx} = ensemble_add_data_struct_row(...
        outdata.data{paths_idx},'subject_id',subid,'session',isess,...
        'run',0,'path_type','spm_outdir','path',spm_outdir);
  end
end
