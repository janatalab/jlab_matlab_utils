function as = repack_formdata(as,ed)
% repackages data from one or more forms into a more useable format
%
% as = repack_formdata(as,ed)
%
% as -- analysis structure (see init_analysis_struct.m for details)
% ed -- experiment data
%
% Multiple instances of a stimulus are associated with unique entries along the 3rd dimension
% of the .data field and the 2nd dimension of the .stimidx field.  Scripts that
% subsequently operate on the as structure can easily locate the multiple
% instances by searching for a specific stimidx in the stimidx field.

% 08/19/05 Petr Janata
% 02/19/06 PJ - Added checking to make sure that there are entries in the form
%               to be processed.
% 03/05/06 PJ - Added code to make sure that the rows in the .num and .txt data
%               matrices correspond to the same subjects.
%               Added writing of datenum information.
% 07/18/06 PJ - Improved handling of stimulus vs non-stimulus numerical
%               responses
% 07/19/06 PJ - Resolved handling of multiple iterations of a stimulus.
% 07/27/06 PJ - Added provisions for two types of timestamps, one that is
%               stimulus based, and one that is question based
% 08/03/06 PJ - Eliminated separate handling of single iteration data since
%               this was causing minor problems.

NUMERIC = 1;

DEBUG = 2;

numeric_types = {'int16','int32','int64','double','enum'};

p = as.params;

nforms = length(as.forms);

iter_cnt_num = [];
iter_cnt_txt = [];

% We need to go through the forms twice.
% The first time is to extract all of the existing subject IDs and then create
% a sorted list of these.  This list will structure the .data tables and will
% ensure that the rows of the .num and .txt data structures are in register.
% It is much easier to ensure that this reasonable assumption is true at this
% stage rather than having to do explicit matching by subject ID in scripts
% that analyze the analysis structure. Note that mismatches in subids may
% persist across elements of an analysis structure, e.g. as(1) vs as(2).

master_sublist = {};
for ifrm = 1:nforms
  cfname = as.forms{ifrm};  % get the current form name
  frm_idx = strmatch(cfname,ed.form_names,'exact');

  % copy the current formdata
  fd = ed.form_data{frm_idx};

  % make sure there are entries to process
  nentries = size(fd.data{1},1);
  if nentries == 0
    continue
  end
  
  % Map our form column constants
  FD = set_form_col_const(fd.vars);

  % Get the list of subjects for whom we have data in this form
  subids = fd.data{FD.SUB_ID};   
  
  % Update the master list
  master_sublist = union(master_sublist,subids);
end  % for ifrm=

% set the subid fields to have the master list
as.num.subid = master_sublist;
as.txt.subid = master_sublist;

for ifrm = 1:nforms
  cfname = as.forms{ifrm};  % get the current form name
  frm_idx = strmatch(cfname,ed.form_names,'exact');
  
  % Make a local copy of the form data
  fd = ed.form_data{frm_idx};
  nentries = size(fd.data{1},1);
  
  if nentries == 0
    continue
  end
  
  % Map our form column constants
  FD = set_form_col_const(fd.vars);

  % Get some basic subject info
  subids = fd.data{FD.SUB_ID};   
  unique_subs = unique(subids);
  nsub = length(unique_subs);
  
  % Determine what the questions are on this form and create appropriate
  % entries in the output structures if they don't yet exist
  
  qid = cat(2,fd.data{[FD.QUEST_ID FD.SUBQUEST_ID]});  % get question IDs
  compqid = qid(:,1)+qid(:,2)/10;

  [unique_quest_ids, quest_idxs] = unique(qid,'rows'); % get the unique questions
  nquest = size(unique_quest_ids,1);
  
  % Create the output question ID field
  out_ids = unique_quest_ids(:,1)+unique_quest_ids(:,2)/10; % qid.sqid
  
  % Retrieve the data formats
  qid_str = sprintf('(question_id=%d AND subquestion=%d) OR ', unique_quest_ids');
  qid_str(end-3:end) = [];  
  mysql_str = sprintf(['SELECT type, data_format_id FROM data_format ' ...
	'RIGHT JOIN question_x_data_format ON' ...
	' data_format.data_format_id=question_x_data_format.answer_format_id ' ...
	'WHERE (%s);'], qid_str);
  [types, dfid] = mysql(p.conn_id,mysql_str);
  
  % Find out which of the entries have to do with a stimulus ID
  stim_mask = ~isnan(fd.data{FD.STIM_ID});
  
  % Determine which ones are numeric
  is_numeric = ismember(types,numeric_types);
  
  for itype = 1:2
    
    % Define a few variables that depend on the data type
    if itype == NUMERIC  % numeric data
      ts = as.num;
      type_mask = is_numeric;
      src_col = FD.RESP_ENUM;  % column in fd.data that we pull data from
      iter_cnt = iter_cnt_num;
    else  % non-numeric data
      ts = as.txt;
      type_mask = ~is_numeric;
      src_col = FD.RESP_TXT;  % column in fd.data that we pull data from
      iter_cnt = iter_cnt_txt;
    end
    
    % Check to see if we need to handle any data of this type
    if ~any(type_mask)
      continue
    end
    
    % Check to see which ones we need to add to the numeric and text output structures
    proc_qid_idxs = find(~ismember(out_ids,ts.qid) & type_mask);
    nproc = length(proc_qid_idxs);
  
    if nproc
      % copy the new question IDs
      insert_idxs = (length(ts.qid)+1):(length(ts.qid)+nproc);
      ts.qid(insert_idxs) = out_ids(proc_qid_idxs);
      
      % copy the new data format IDs (for decoding enums)
      if itype == 1
	ts.dfid(insert_idxs) = dfid(proc_qid_idxs);
      end
      
      % copy the question text
      tmp = cat(2,fd.data{[FD.QUEST_TXT FD.SUBQUEST_TXT]});
      for iproc = 1:nproc
	if ~strcmp(tmp{quest_idxs(iproc),1},tmp{quest_idxs(proc_qid_idxs(iproc)),2})
	  ts.qtxt{insert_idxs(iproc)} = ...
	      cell2str(tmp(quest_idxs(proc_qid_idxs(iproc)),:),'\n');
	else
	  ts.qtxt{insert_idxs(iproc)} = tmp{quest_idxs(proc_qid_idxs(iproc)),1};
	end
      end
      
      % Initialize some iteration handling variables
%       if isempty(ts.niter)
% 	last_idx = nsub;
%       else
% 	last_idx = size(ts.niter,1);
%       end
%       ts.niter(last_idx,insert_idxs) = 0  % Total number of iterations/subject/question
%       iter_cnt(last_idx,insert_idxs) = 0;
    end % if nproc
    
    %
    % COPY THE DATA FOR EACH SUBJECT.
    %
    % Exactly how this happens depends on the how we want to handle multiple
    % responses to the same questions, as happens when the same form is
    % presented multiple times. This can happen either in the context of
    % non-stimulus based forms, e.g. PANAS, or with forms that are associated
    % with a stimulus ID.
    %
    % The easiest way to handle this might be to place repetitions into a 3rd
    % dimension in the data field.
    %
    
    for isub = 1:nsub
      sid = unique_subs{isub};
      
      % Check to see if we already have an entry for this subject and create if
      % necessary.
      % This conditional code should not be entered given the changes made of
      % initializing the subid lists in the .num and .txt structures to be the
      % same (same elements, same order).
      
      row_idx = find(strcmp(ts.subid,sid));
      if isempty(row_idx)
	if DEBUG == 2
	  fprintf('Creating Type %d entry for subject: %s\n', itype, sid);
	end
	ts.subid{end+1} = sid;
	row_idx = length(ts.subid);
      end

      % Make sure our iteration counters are up to date
      if (size(iter_cnt,1) < row_idx) | (size(iter_cnt,2) < length(ts.qid))
	ts.niter(row_idx,length(ts.qid)) = 0;
	iter_cnt(row_idx,length(ts.qid)) = 0;
      end
      
      % Determine which entries in the form data belong to this subject
      submask = strcmp(subids,sid);
      
      % Determine destination columns for the data. compqid(submask) pulls out
      % the question IDs that we have responses for for this subject, and
      % ts.qid is a vector that maps each question ID to a column in the output
      % data matrix.  have_dest indicates which of the responses can be copied,
      % and dest_idx gives the appropriate column in the output data matrix.

      destmask = zeros(size(submask));
      [have_dest,dest_idx] = ismember(compqid(submask),ts.qid);

      % Make sure we have some data to enter
      if any(have_dest)
	destmask(submask) = have_dest;
	
	if itype == NUMERIC & (DEBUG == NUMERIC)
	  [dest_idx(have_dest) fd.data{src_col}(submask&destmask)]
	  [sum(submask&destmask) length(dest_idx(have_dest))]
	end
	
	% update the iteration count
	dest_cols = unique(dest_idx(have_dest));
	ndest = length(dest_cols);
	if ndest > 1
	  ts.niter(row_idx,dest_cols) = ...
	      ts.niter(row_idx,dest_cols) + ...
	      hist(dest_idx(have_dest),dest_cols);
	else
	  ts.niter(row_idx,dest_cols) = ...
	      ts.niter(row_idx,dest_cols) + ...
	      sum(have_dest);
	end
	
	% Now copy the data.
	
	  % Now loop over destination columns (individual questions)
	  for idest = 1:ndest
	    col_idx = dest_cols(idest);
	    % Find which responses pertain to this question.
	    colmask = compqid == ts.qid(col_idx);
	    
	    curr_mask = colmask&submask;
	    clear colmask
	    
	    % check to see if we are dealing with stimulus IDs
	    is_stim = any(fd.data{FD.STIM_ID}(curr_mask&stim_mask));
	    
	    if is_stim
	      curr_mask = curr_mask&stim_mask;
	    end
	  
	    % Figure out how many responses we have to this particular question
	    nresp = sum(curr_mask);
	    
	    % Copy the data. If we aren't dealing with a stimulus, then it is easy
	    % because we don't have to go into specific stim_id slots, and can
	    % just pack things into the last place we left off.  
	    % 
	    % NOTE: This is a somewhat weak solution, meaning that there might
	    % be situations where we want things aligned by iteration even if we
	    % aren't dealing with a stimulus.  So, the handling of this may need
	    % to be changed in the future.
	    if ~is_stim
	      % Figure out where the responses are going to go in the 3D data
	      % field.
	      rep_idxs = ...
		  iter_cnt(row_idx,col_idx)+1:iter_cnt(row_idx,col_idx)+nresp;
	      
	      % update the iteration count
	      iter_cnt(row_idx,col_idx) = max(rep_idxs);
	      
	      % copy the data
	      ts.data(row_idx,col_idx,rep_idxs) = ...
		  fd.data{src_col}(curr_mask);
	      
	      ts.datenum.by_question(row_idx,col_idx,rep_idxs) = fd.data{FD.DATE_TIME}(curr_mask);
	      %diff(fd.data{FD.RESP_ID}(curr_mask)) < 0
	    else
	      % Handle each stimulus individually.  
	      % We have to go through some extra hoops to preserve the temporal
	      % order in which they occurred.
	      % 
	      % HANDLING OF MULTIPLE ITERATIONS OF A STIMULUS
	      %
	      % Multiple iterations of a stimulus ID are handled by simply
	      % adding another instance of the stimulus ID to the .data and
	      % .stimidx matrices (the 3rd and 2nd dimensions, respectively),
	      % rather than creating an extra dimension to handle multiple
	      % stimulus iterations.  The advantage of this scheme is that it
	      % preserves the temporal ordering of stimulus presentation in
	      % these variables. Thus, unique instances of a stimulus are now
	      % identified by a combination of their stimulus ID and timestamp.
	      
	      [curr_stim_ids,idx1,idx2] = ...
		  unique(fd.data{FD.STIM_ID}(curr_mask));
	      curr_mask_idxs = find(curr_mask);
	      
	      % Make sure we have no NaNs among the stim ids
	      bad_ids = isnan(curr_stim_ids);
	      if any(bad_ids)
		warning(sprintf('Found %d NaNs among the stim IDs. subid=%s\n', sum(bad_ids), subid))
	      end
	      
	      % sort them into original temporal order
	      curr_stim_ids = curr_stim_ids(idx2);
	      curr_stim_times = fd.data{FD.DATE_TIME}(curr_mask_idxs);
	      
	      % Make sure that there are no duplicate stimulus submissions
	      dup_idxs = find((diff(curr_stim_ids)==0) & ...
		  (diff(curr_stim_times)==0))+1;

	      if ~isempty(dup_idxs)
		warning(sprintf(['\n%d duplicated responses with same stimulus ' ...
		      'ID (%s) and timestamp: subject (%s)\n'], length(dup_idxs), sprintf('%d,', curr_stim_ids(dup_idxs)), sid))
		curr_stim_ids(dup_idxs) = [];
		curr_stim_times(dup_idxs) = [];
		nresp = nresp-1;
	      end
	      
	      num_presented_stims = length(curr_stim_ids);
	      
	      % issue a warning if there are more responses than stimuli,
	      % i.e. some stimulus was responded to more than once
	      if num_presented_stims ~= nresp
		warning(sprintf('Encountered %d stimuli and %d responses\n', num_presented_stims, nresp))
	      end
	      
	      for istim = 1:num_presented_stims
		curr_stim_id = curr_stim_ids(istim);
		
		% Check to see if this particular stim_id already exists in the
		% master stimulus list. If not, add it to the list.
		master_stim_idx = find(as.stims.ids==curr_stim_id);
		
		if isempty(master_stim_idx)
		  master_stim_idx = length(as.stims.ids)+1;
		  
		  % Copy the stim ID to this location
		  as.stims.ids(master_stim_idx) = curr_stim_id;
		end % if isempty(master_stim_idx)
		
		% Now, check to see if this particular master_stim_idx already
		% exists in the subject's stimulus list, and if it doesn't,
		% pack it into the first available slot
		
		% Make sure we've started a row in the stimidx matrix for this subject
		if size(ts.stimidx,1) < row_idx
		  ts.stimidx(row_idx,1) = 0;
		end
		
		ts.stimidx(row_idx,istim) = master_stim_idx;
		ts.datenum.by_stim(row_idx,istim) = curr_stim_times(istim);
		
		% make a mask for this instance of the stimulus
		curr_stimmask = (fd.data{FD.STIM_ID} == curr_stim_id) & ...
		    (fd.data{FD.DATE_TIME} == curr_stim_times(istim)); 
		
		% Copy the data to this location
		tmpdata = fd.data{src_col}(curr_mask&curr_stimmask);
		
		% Only 1 data value should be returned, though in the case of a
		% doubly submitted response, there will be more than
		% one. Therefore, explicitly take the first.
		ts.data(row_idx,col_idx,istim) = tmpdata(1);
		
		% update the iteration count
		iter_cnt(row_idx,col_idx) = iter_cnt(row_idx,col_idx)+1;
	      end % for istim=
	    end % if ~is_stim
	  end % for idest=

      end % if any(have_dest)
    end % for isub=

    % put the temporary structure back into the proper output structure
    if itype == NUMERIC
      as.num = ts; 
      iter_cnt_num = iter_cnt;
    else 
      as.txt = ts; 
      iter_cnt_txt;
    end
    
  end % for itype=
end % for ifrm=
