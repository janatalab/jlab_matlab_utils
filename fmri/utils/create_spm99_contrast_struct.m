function [contrasts] = create_contrast_struct(sinfo,cont_name,cont_type,incl_code)
% [contrasts] = create_contrast_struct()
%
% Utility function for populating the contrasts structure that is used by SPM99 
%
% Previously, this script was part of the *_contrasts.m file associated with an experiment

% 12/31/03 PJ
% 09/09/04 PJ  Added analysis path variable

global rootpath model_proto_name proc_subs analysis_path

subject_ids = cellstr(char(sinfo(:).id));
nsub_proc = length(proc_subs);
nsub = length(sinfo);
nc = length(cont_name);

for isub = 1:nsub_proc
  co = {};
  sub_idx = proc_subs(isub);
  
  subpath = fullfile(rootpath, subject_ids{sub_idx}, analysis_path, model_proto_name);  
  
  %
  % Get the column names for this subject
  %
  % Selection of column names can happen in one of two different ways.
  % If conditions are modeled with several basis functions, and the contrast is
  % specified by the names of the conditions, then the column
  % indices must be grabbed from the Sess structure array.
  %
  % If contrasts are specified with individual column names, the column indices
  % are retrieved by matching the names specified in the contrast with the
  % names in xX.Xnames
  %
  % Those names ending with a number in parentheses refer to individual
  % columns/basis functions, as listed in xX.Xnames
  %
  % Column indices associated with names that don't end with a number in
  % parentheses are *safely* retrieved from the Sess structure. Note that there
  % is a Sess structure for each run
  
  x = load(fullfile(subpath,'SPM.mat'));
  cond_names_detail = x.xX.Xnames;
  Sess = x.Sess;
  clear x
  
  % Deal with the Sess mess
  cond_names_grouped = [];
  tot_ngroup = 0;
  for isess = 1:length(Sess)
    cond_names_grouped = [cond_names_grouped Sess{isess}.name];
    ngroup = length(Sess{isess}.name);
    for igroup = 1:ngroup
      tot_ngroup = tot_ngroup+1;
      grouping_idxs{tot_ngroup} = Sess{isess}.col(Sess{isess}.ind{igroup});
    end
  end
  
  cond_names_detail = strvcat(cond_names_detail);
  cond_names_detail = cond_names_detail(:,7:end);	% strip leading session identifiers
  cond_names_detail = cellstr(cond_names_detail);

  good_cont_vect = [];
  ngood = 0;
  for icont = 1:nc
    switch cont_type{icont}
      case 'T'
	cond_a_idx = [];
	cond_b_idx = [];
	% First, try to grab the column index from the detailed list
	% This will pick up any regressors and condition components
	cond_a_idx = find(ismember(cond_names_detail, incl_code{icont,1}));
	cond_b_idx = find(ismember(cond_names_detail, incl_code{icont,2}));
	
	% Now, grab column indices based on group specifiers
	tmp_idx = cat(2,grouping_idxs{find(ismember(cond_names_grouped, ...
	      incl_code{icont,1}))});
	cond_a_idx = [cond_a_idx tmp_idx];
	
	tmp_idx = cat(2,grouping_idxs{find(ismember(cond_names_grouped, ...
	      incl_code{icont,2}))});
	cond_b_idx = [cond_b_idx tmp_idx];
	
	cond_a_vect = zeros(1, length(cond_names_detail));
	cond_b_vect = zeros(1, length(cond_names_detail));
	
	cond_a_vect(cond_a_idx) = 1;
	cond_b_vect(cond_b_idx) = 1;
	
	% Weight the contrasts appropriately, i.e. have sum=0
	suma = sum(cond_a_vect);
	sumb = sum(cond_b_vect);
	
	if suma & sumb
	  ngood = ngood+1;
	  good_cont_vect(icont) = 1;
	  
	  weight_b = suma/sumb*-1;
	  
	  cond_b_vect = cond_b_vect * weight_b;
	  
	  co{ngood} = cond_a_vect + cond_b_vect;
	elseif suma
	  %disp('Only numerator was specified in T contrast')
	  ngood = ngood+1;
	  good_cont_vect(icont) = 1;
	  co{ngood} = cond_a_vect;
	else
	  good_cont_vect(icont) = 0;
	  disp(sprintf('Skipping T contrast (%s): num_num: %d; num_denom: %d', cont_name{icont},suma, sumb))
	end
	
      case 'F'
	col_idx = [];
	
	col_idx = [col_idx find(ismember(cond_names_detail, incl_code{icont,1}))];
	
	tmp_idx = cat(2,grouping_idxs{find(ismember(cond_names_grouped, ...
	      incl_code{icont,1}))});
	
	col_idx = [col_idx tmp_idx];
	
	cols = zeros(1,length(cond_names_detail));
	cols(col_idx) = 1;
	
	if any(cols)
	  ngood = ngood+1;
	  good_cont_vect(icont) = 1;
	  co{ngood} = full(sparse(1:sum(cols),find(cols),1));
	  co{ngood}(end,length(cond_names_detail)) = 0;
	else
	  good_cont_vect(icont) = 0;
	  disp(sprintf('Skipping F contrast (%s)', cont_name{icont}))
	end
	
      otherwise
	error(sprintf('Unknown contrast type: %s', cont_type{icont}))
	
    end % switch cont_type{icont}
    
  end % for icont

  good_idxs = find(good_cont_vect);
  contrasts(sub_idx).names = cont_name(good_idxs);
  contrasts(sub_idx).types = cont_type(good_idxs);
  contrasts(sub_idx).values = co;
end
