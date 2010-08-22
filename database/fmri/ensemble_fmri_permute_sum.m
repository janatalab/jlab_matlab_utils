function outdata = ensemble_fmri_permute_sum(indata,params)

% generate summary image from permuted fmri models
% 
%   outdata = ensemble_fmri_permute_sum(indata,params)
% 
% This function loads all permutation probability images provided to it
% (from the 'permprob' output of ensemble_fmri_build_model_l1), applies a
% probability cutoff to each image (params.permute_sum.thresh), then masks
% out all voxels not belonging to a cluster of sufficient size
% (params.permute_sum.k). It then sums across all subject images, saves
% that summary image to disk, and outputs the path of that image.
% 
% You can use ensemble_concat_datastruct to concatenate the permprob
% datasets from your subjects into one datastruct, and pass that datastruct
% to this function. You are free to include datastructs from multiple
% models (just make sure you include all subjects that you are interested
% in for a given model). In the case that you provide data from multiple
% models, this function will generate separate summary images for separate
% models.
% 
% REQUIRES
%   indata
%       permprob - permutation probability maps output by
%           ensemble_fmri_build_model_l1.m
%   params.permute_sum
%       .thresh - permutation probability inclusion threshold. This
%           is the probability of a voxel containing real signal, therefore
%           is probably a value of 0.95 or above (as opposed to 0.05 or
%           below). default: 0.95
%       .k - cluster size threshold. default: 40
% 
% RETURNS
%   outdata - ensemble data struct, type "permute_sum", containing the model
%       id, path, and the number of images that were used to create each
%       permutation sum image
% 
% FB 2010.08.20

global r;

outdata = '';

%% return output directory?
if (iscell(indata) && ~isempty(indata) && isfield(indata{1},'task') && ...
        ~isempty(strmatch('return_outdir',indata{1}.task))) || ...
        (isstruct(indata) && isfield(indata,'task') && ...
        ~isempty(strmatch('return_outdir',indata.task)))
  if isfield(params,'paths') && isfield(params.paths,'analpath')
    outdata = params.paths.analpath;
  end
  return
end

%% init vars
outdata = ensemble_init_data_struct();
outdata.type = 'permute_sum';
outdata.vars = {'model_id','path','nimg'};
ocol = set_var_col_const(outdata.vars);
outdata.data{ocol.model_id} = [];
outdata.data{ocol.path} = {};
outdata.data{ocol.nimg} = [];

r = init_results_struct;

r.type = 'fmri_anal';  % Identify the type of this reporting instance
r.report_on_fly = 1;

% Parse out the input data
for idata = 1:length(indata)
  if isfield(indata{idata},'type')
      switch indata{idata}.type
        case 'permprob'
          permdata = indata{idata};
          permcol = set_var_col_const(permdata.vars);
          modids = unique(permdata.data{permcol.model_id});
          nmod = length(modids);
      end
  end
end

% check for required vars, quit if they can't be found
check_vars = {'permdata'};
check_required_vars;

% get perm sum parameters
try t = params.permute_sum.thresh; catch t = 0.95; end % probability threshold
try k = params.permute_sum.k; catch k = 40; end        % cluster size threshold

%% generate permutation analysis summary image
for k=1:nmod
    
  % init vars
  mid = modids(k);
  mfilt.include.all.model_id = mid;
  mdata = ensemble_filter(permdata,mfilt);
  
  nimg = length(mdata.data{1});
  
  V = spm_vol(mdata.data{permcol.path}{1});
  Z = zeros(V.dim);

  % load images, filter, sum clusters
  for l=1:nimg
    V = spm_vol(mdata.data{permcol.path}{l});
    Y = spm_read_vols(V);
    Y = Y > t;
    Y = fmri_clust_filt(Y,k);
    
    Z = Z + Y;
  end % for l=1:nsub

  % save img to file
  opath = fullfile(params.paths.analpath,'spm','group',...
      sprintf('model_%d',mid));
  check_dir(opath,0,1);
  Vout = V;
  Vout.fname = fullfile(opath,...
      sprintf('mod%d_permute_sum_%dimgs.img',mid,nimg));
  spm_write_vol(Vout,Z);

  % save to outdata
  outdata = ensemble_add_data_struct_row(outdata,...
      'model_id',mid,'path',Vout.fname,'nimg',nimg);
end % for k=1:nmod

fprintf(1,'DONE\n');
