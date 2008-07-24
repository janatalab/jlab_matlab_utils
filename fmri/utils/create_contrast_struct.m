function [consess] = create_contrast_struct(spmmat_fname,cinfo)
% [consess] = create_contrast_struct()
%
% Utility function for populating the consess field of a job structure that is
% used by the SPM5 contrast manager
%
% Name of SPM.mat file to which the .Xcon structure will ultimately be added
%
% cinfo is a structure with fields
%   .cont_name - teh contrast name
%   .cont_type - the contrast type (T,F)
%   .incl_code - the conditions to include in the contrast
%
% Multiple contrasts are specified as multiple elements in a cinfo structure array

% 02/23/06  Petr Janata
% 08/10/06  PJ - adapted from earlier version. Made cinfo to be a structure
%                array rather than a single structure filled with cell arrays

if ~exist(spmmat_fname,'file')
  error(sprintf('Could not locate SPM.mat file: %s',spmmat_fname))
end

%
% Load information from the SPM.mat file about the regressors in the design matrix
% 
load(spmmat_fname);
cond_names_detail = SPM.xX.name;
%Sess = SPM.Sess;
clear SPM

ncond = length(cond_names_detail);
cond_names_stripped = cell(1,ncond);
% Strip out the basic condition names that we want to match against
for icond = 1:ncond
  condstr = cond_names_detail{icond};
  start_idx = min(find(isspace(condstr)))+1;
  if isempty(start_idx)
    start_idx = 1;
  end
  stop_idx = min(regexp(condstr,'[\^*]'))-1; % Search for a caret and asterisk
  if isempty(stop_idx)
    stop_idx = length(condstr);
  end
  cond_names_stripped{icond} = condstr(start_idx:stop_idx);
end

% Initialize some variables
good_cont_vect = [];
ngood = 0;
nc = length(cinfo); % Number of contrasts

tcon = struct('name','','convec',[]);
fcon = struct('name','','convec',[]);
consess = {};

% Loop over the list of contrasts
for ic = 1:nc
  cont_type = cinfo(ic).cont_type;
  cont_name = cinfo(ic).cont_name;
  incl_code = cinfo(ic).incl_code;
  switch cont_type
    case 'T'
      cond_a_idx = [];
      cond_b_idx = [];
      % First, try to grab the column index from the detailed list
      % This will pick up any regressors and condition components
      cond_a_idx = find(ismember(cond_names_stripped, incl_code{1}));
      cond_b_idx = find(ismember(cond_names_stripped, incl_code{2}));
    
      cond_a_vect = zeros(1, length(cond_names_stripped));
      cond_b_vect = zeros(1, length(cond_names_stripped));
      
      cond_a_vect(cond_a_idx) = 1;
      cond_b_vect(cond_b_idx) = 1;
      
      % Weight the contrasts appropriately, i.e. have sum=0
      suma = sum(cond_a_vect);
      sumb = sum(cond_b_vect);
      
      if suma & sumb
	ngood = ngood+1;
	weight_b = suma/sumb*-1;
	
	cond_b_vect = cond_b_vect * weight_b;
	
	convec = cond_a_vect + cond_b_vect;
      elseif suma
	%disp('Only numerator was specified in T contrast')
	ngood = ngood+1;
	convec = cond_a_vect;
      elseif sumb
	ngood = ngood+1;
	convec = cond_b_vect*-1;
      else
	good_cont_vect(ic) = 0;
	disp(sprintf('Skipping T contrast (%s): num_num: %d; num_denom: %d', cont_name,suma, sumb))
	continue
      end
      consess{ngood}.tcon.name = cont_name;
      consess{ngood}.tcon.convec = convec;
      
    case 'F'
      col_idx = find(ismember(cond_names_stripped, incl_code{1}));
      
      if 0
      tmp_idx = cat(2,grouping_idxs{find(ismember(cond_names_grouped, ...
	    incl_code{1}))});
      col_idx = [col_idx tmp_idx];
    end
      
      cols = zeros(1,length(cond_names_stripped));
      cols(col_idx) = 1;
      
      if any(cols)
	ngood = ngood+1;
	convec = full(sparse(1:sum(cols),find(cols),1));  % actually a matrix
	convec(end,length(cond_names_stripped)) = 0; % pad out with zeros
	consess{ngood}.fcon.name = cont_name;
	for irow = 1:size(convec,1)
	  consess{ngood}.fcon.convec{irow} = num2str(convec(irow,:));
	end
      else
	good_cont_vect(ic) = 0;
	disp(sprintf('Skipping F contrast (%s)', cont_name))
	continue
      end
      
    otherwise
      error(sprintf('Unknown contrast type: %s', cont_type))
      
  end % switch cont_type{ic}
  
end % for ic

return
