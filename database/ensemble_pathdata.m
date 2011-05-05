function outdata = ensemble_pathdata(indata,defs)

% returns an ensemble pathdata struct, given defs.paths and defs.sinfo
% 
%   outdata = ensemble_pathdata({},defs)
% 
% REQUIRES
%   indata - empty cell array
%   defs.paths - a structure containing path data
%   defs.sinfo - a structure containing subject information
% 
% FB 2011.04.29

%% init outdata
outdata = ensemble_init_data_struct();
outdata.type = 'paths';
outdata.vars = {'subject_id','session','run','path_type','path'};
pathcol = set_var_col_const(outdata.vars);
outdata.data{pathcol.subject_id} = {};
outdata.data{pathcol.session} = {};
outdata.data{pathcol.run} = [];
outdata.data{pathcol.path_type} = {};
outdata.data{pathcol.path} = {};

%% iterate over subjects
nsub = length(defs.sinfo);
for k=1:nsub
  subid = defs.sinfo(k).id;
  
  % Deal with directory infrastructure
  sub_indir = fullfile(defs.paths.inroot, subid);
  sub_outdir = fullfile(defs.paths.outroot, subid);
  check_dir(sub_outdir,1);
	
  % save paths to the path output data structure
  % subject source file directory
  outdata = ensemble_add_data_struct_row(outdata,'subject_id',subid,...
      'session',{''},'run',0,'path_type','sub_indir','path',sub_indir);
	
  % subject destination file directory (often the same as the source)
  outdata = ensemble_add_data_struct_row(outdata,'subject_id',subid,...
      'session',{''},'run',0,'path_type','sub_outdir','path',sub_outdir);

  nsess = length(defs.sinfo(k).sessinfo);
  for is=1:nsess
    sessid = defs.sinfo(k).sessinfo(is).id;
    
    % define paths
    sess_indir = fullfile(sub_indir,sessid);
    sess_outdir = fullfile(sub_outdir,sessid);
    behav_indir = fullfile(sess_indir,'behavior');
    behav_outdir = fullfile(sess_outdir,'behavior');
    physio_indir = fullfile(sess_indir,'physio');
    physio_outdir = fullfile(sess_outdir,'physio');
    
    % check paths
    check_dir(sess_indir,1);
    check_dir(sess_outdir,1);
    check_dir(behav_indir,1);
    check_dir(behav_outdir,1);
    check_dir(physio_indir,1);
    check_dir(physio_outdir,1);
  
    % save paths
    outdata = ensemble_add_data_struct_row(outdata,'subject_id',subid,...
        'session',sessid,'run',0,'path_type','sess_indir',...
        'path',sess_indir);
    outdata = ensemble_add_data_struct_row(outdata,'subject_id',subid,...
        'session',sessid,'run',0,'path_type','sess_outdir',...
        'path',sess_outdir);
    outdata = ensemble_add_data_struct_row(outdata,'subject_id',subid,...
        'session',sessid,'run',0,'path_type','behav_indir',...
        'path',behav_indir);
    outdata = ensemble_add_data_struct_row(outdata,'subject_id',subid,...
        'session',sessid,'run',0,'path_type','behav_outdir',...
        'path',behav_outdir);
    outdata = ensemble_add_data_struct_row(outdata,'subject_id',subid,...
        'session',sessid,'run',0,'path_type','physio_indir',...
        'path',physio_indir);
    outdata = ensemble_add_data_struct_row(outdata,'subject_id',subid,...
        'session',sessid,'run',0,'path_type','physio_outdir',...
        'path',physio_outdir);
  end % for is=1:nsess
end % for k=1:nsub
