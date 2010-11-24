function outdata = ensemble_fmri_load_epi(indata,defs)
% Loads EPI files from a directory
%
% outdata = ensemble_fmri_load_epi(indata,defs);
%
% The information needed to construct the path to the files is specified in
% defs:
%
% defs.paths - general path information about the project
% defs.sinfo - subject information, potentially containing information
%              about multiple sessions.
% defs.sessinfo - the session information struct.  The sessinfo.id field 
%               must correspond to a sub-directory within
%               the subject directory.  If sessinfo is not specified or
%               is left empty, all sessions specified in sinfo.sessinfo
%               will be used in accordance with their
%               sinfo.sessinfo.use_session flag.  The default is to not use
%               a session.
%
% defs.load_runs - runs within a session to load. If this is not specified,
%               sessinfo.use_runs will be used.
%
% defs.file_prefix - the search string used to specify the files to load,
%               e.g. 'sw' to specify smoothed, normalized files
%
% outdata - an EPI data structure is returned

% 18Nov2010 Petr Janata

sinfo = defs.sinfo;
subid = sinfo.id;

if 0  % nested version
outdata = ensemble_init_data_struct();
outdata.type = 'epi';

outdata.vars{end+1} = 'epi';
epi_idx = length(outdata.vars);

outdata.data{epi_idx} = ensemble_init_data_struct();
outdata.data{epi_idx}.type='epi';
outdata.data{epi_idx}.vars = {'subject_id','session',...
	'ensemble_id','run','path'};
epiCols = set_var_col_const(outdata.data{epi_idx}.vars);

outdata.data{epi_idx}.data{epiCols.subject_id} = {};
outdata.data{epi_idx}.data{epiCols.session} = [];
outdata.data{epi_idx}.data{epiCols.ensemble_id} = [];
outdata.data{epi_idx}.data{epiCols.run} = [];
outdata.data{epi_idx}.data{epiCols.path} = {};
end

% top-level version
outdata = ensemble_init_data_struct();

outdata.type='epi';
outdata.vars = {'subject_id','session',...
	'ensemble_id','run','path'};
epiCols = set_var_col_const(outdata.vars);

outdata.data{epiCols.subject_id} = {};
outdata.data{epiCols.session} = [];
outdata.data{epiCols.ensemble_id} = [];
outdata.data{epiCols.run} = [];
outdata.data{epiCols.path} = {};

if isfield(defs, 'sessid')
	sessidx = find(ismember({sinfo.sessinfo(:).id}, defs.sessid));
	if isempty(sessidx)
		error('Failed to find session (%s) for subject (%s)', defs.sessid, subid)
	end
else
	sessidx = 1:length(sinfo.sessinfo);
end

nsess = length(sessidx);
for isess = 1:nsess
	idx = sessidx(isess);
	sesspath = fullfile(defs.paths.outroot, subid, sinfo.sessinfo(idx).id);
	
	try load_runs = defs.load_runs; catch load_runs = []; end
	if isempty(load_runs)
		load_runs = sinfo.sessinfo(idx).use_runs;
	end
	
	for irun = 1:length(load_runs)
		runidx = load_runs(irun);
		runpath = fullfile(sesspath, 'epi', sprintf('run%d', runidx));
		
		% get the filelist
		flist = get_spm_flist(runpath, defs.file_prefix);
		
		% populate the output structure
		for ifl=1:size(flist,1)
			outdata = ensemble_add_data_struct_row(...
				outdata,'subject_id',subid,...
				'session',idx,...
				'ensemble_id', sinfo.sessinfo(idx).ensemble_id,...
				'run',runidx,...
				'path',flist(ifl,:));
		end
	end
end