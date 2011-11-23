function outvar = ensemble_form_table(data_st,params)
% Prints responses to questions on a form out in tabular format
%
% ensemble_form_table(data_st,params)
%
% 05/01/07 Petr Janata
% 23Nov2011 PJ - improved conn_id handling, fixed handling of array of
% forms
outvar = [];

try conn_id = params.conn_id; catch
	if isfield(params,'mysql')
		conn_id = params.mysql.conn_id;
	else
		try conn_id = params.ensemble.conn_id; catch
			conn_id = [];
		end
	end
end

% Make sure we have compqids
data_st = ensemble_check_compqid(data_st);

% Set column indexing variable
cols = set_var_col_const(data_st.vars);

form_ids = unique(data_st.data{cols.form_id});
nforms = length(form_ids);

% Extract form metadata
form_data = mysql_extract_metadata('table','form', ...
    'form_id',form_ids, ...
    'conn_id', conn_id);

dataformat = mysql_extract_metadata('table','data_format','conn_id',conn_id);

% Generate subject masks
[sub_mask_mtx, subids] = make_mask_mtx(data_st.data{cols.subject_id});
nsub = length(subids);

for iform = 1:nforms
  form_id = form_data(iform).form_id;
  form_mask = data_st.data{cols.form_id} == form_id;
  
  % Format things so that we can call ensemble_display_table
  table_data = ensemble_init_data_struct;
  table_data.vars = {'data','column_labels','column_formats'};
  tblcols = set_var_col_const(table_data.vars);

  for isub = 1:nsub
    
    nquest = length(form_data(iform).question);
    for iq = 1:nquest
      q = form_data(iform).question(iq);
      compqid = q.compqid;
      qid_mask = data_st.data{cols.compqid} == compqid;
      
      table_data.data{tblcols.column_labels}{iq} = q.question_text;
      
      % Make a composite mask
      compmask = form_mask & qid_mask & sub_mask_mtx(:,isub);
      if ~any(compmask)
				continue
      end
      
      % Determine what type of data the question is
      df_idx = find([dataformat.data_format_id]==q.data_format_id);
      qtype = dataformat(df_idx).type;
      switch qtype
	case 'enum'
	  qvals = dataformat(df_idx).enum_values;
	  qdata = data_st.data{cols.response_enum}(compmask);
    switch q.html_field_type
      case 'checkbox'
        qdata = logical(data2bitmask(qdata));
      otherwise
        qdata = enum2data(qdata);
    end
	  table_data.data{tblcols.data}{iq}{isub} = qvals(qdata); 
	case {'varchar','text'}
	  table_data.data{tblcols.data}{iq}(isub) = ...
	      data_st.data{cols.response_text}(compmask);
	case {'int16','int32','int64','double'}
	  table_data.data{tblcols.data}{iq}(isub) = ...
	      str2double(data_st.data{cols.response_text}(compmask));
	otherwise
	  table_data.data{tblcols.data}{iq}{isub} = ...
	      data_st.data{cols.response_text}(compmask);
      end
      
      % Set the column format
      switch qtype
	case {'enum','varchar','text'}
	  table_data.data{tblcols.column_formats}{iq} = '%s';
	case {'int16','int32','int64'}
	  table_data.data{tblcols.column_formats}{iq} = '%d';
	otherwise
	  table_data.data{tblcols.column_formats}{iq} = '%1.4f';
      end
    end % for iquest =
    
  end % for isub=
  
  % Call the table display function
  ensemble_display_table(table_data, params);
  
end % for iform=
