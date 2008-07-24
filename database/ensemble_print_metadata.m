function ensemble_print_metadata(data_str,params)
% Displays metadata obtained from mysql_extract_metadata
%
%
% ensemble_print_metadata(data_str,params)
%
% data_str is a data structure returned from mysql_extract_metadata
% either at the experiment, form, or question level
%
% params is a structure with display and tables subfields
% params may be left out, in which case default values will be used
% possible subfields of params.tables.report are write2file, fname, and columns
%
%    (e.g. params.report.tables.write2file = 1;
%          params.report.tables.fname      = '/tmp/mymeta.csv';
%          params.report.tables.columns =
%          {{'experiment_id','experiment_title'},...
%           {'form_id','form_name'},...
%           {'question_id','subquestion','question_text','enum_values'}};
%          
% Each cell of the 'columns' field is a cell array of strings, which
% correspond to the fields of each table to print to screen
% or file. They should all be fieldnames of the corresponding
% table(s) with one exception: the string 'answer_format' can be
% used instead of 'enum_values'. Either value produces the same
% output: a column with an "Answer Format" heading that either has
% the enum values for a particular question, or the answer type
% (e.g. text,integer,double, etc.).




DEFAULT_EXPERIMENT_META = {'experiment_id','experiment_title','response_table'};
DEFAULT_FORM_META = {'form_id','form_name','header'};
DEFAULT_QUESTION_META = {'question_id','subquestion','question_text','answer_format'};

%check if write2file parameter was set. If not, then print to screen;
try write2file = params.report.tables.write2file; catch write2file = 0, end
%if write2file was set, open the file and set the fid  
if write2file
  try fid = params.report.tables.fid;
  catch
    fid = fopen(params.report.tables.fname,'wt');
    params.report.tables.fid = fid;
    if fid == -1
      error(sprintf('Problem opening logfile: %s\n',params.report.tables.fname))
    end
    fprintf('Writing tables to file: %s\n', ...
	    params.report.tables.fname);
    closeFid = 1;
  end
else
  fid = 1;
  params.report.tables.fid = fid;
end

%check if column names were set. If not, set to default field names
if(~isfield(params.report.tables,'columns') | isempty(params.report.tables.columns))
  if(isfield(data_str(1),'experiment_id'))
    params.report.tables.columns = {DEFAULT_EXPERIMENT_META,DEFAULT_FORM_META, DEFAULT_QUESTION_META};
  elseif(isfield(data_str(1),'form_id'))
    params.report.tables.columns = {DEFAULT_FORM_META, DEFAULT_QUESTION_META};
  elseif(isfield(data_str(1),'question_id'))
    params.report.tables.columns = DEFAULT_QUESTION_META;
  end
end

%check for table type by the primary key used. There is no table
%name in the meta structure, so currently, this is the only way to
%verify the table type. Printout format depends on table type.
if(isfield(data_str(1),'experiment_id'))
  nForms = length(data_str.form);
  fprintf(fid,'EXPERIMENT TITLE "%s"\n',data_str.experiment_title);
  fprintf(fid,'%d FORMS\n\n',nForms);
  %call the subfunction to print the table
  ensemble_print_table(data_str,params);
  %set the childParams structure, which is the same except the
  %column names should be in the first element of the cell array
  childParams = params;
  childParams.report.tables.columns = params.report.tables.columns(2:end);

  for iForm = 1:nForms
    ensemble_print_metadata(data_str.form(iForm),childParams);
  end
  
elseif(isfield(data_str(1),'form_id'))

  nQuestions = length(data_str.question);
  if(nQuestions == 1)
    qString = 'question';
  else
    qString = 'questions';
  end

  fprintf(fid,'\nForm Name "%s"\n',data_str.form_name);
  fprintf(fid,'%d %s\n\n',nQuestions,qString);
  ensemble_print_table(data_str,params);
  
  if(nQuestions > 0)
    fprintf(fid,'\n');
    childParams = params;
    childParams.report.tables.columns = params.report.tables.columns(2:end);
    ensemble_print_metadata(data_str.question,childParams);
  end
  
elseif(isfield(data_str(1),'question_id'))
  ensemble_print_table(data_str,params);
end

if(write2file & exist('closeFid','var'))
  fclose(fid);
end

return




function ensemble_print_table(data_st,params)
%

fid = params.report.tables.fid;

%ASCII codes for LF, CR, and TAB
%any occurrence of these need to be removed from string
%fields. If these occur in any field in the database, they 
%are only a product of typing style during data entry, since 
%Ensemble only prints carriage returns and tabs in HTML format.
LF_ASCII = 10;
CR_ASCII = 13;
TAB_ASCII = 9;

%these variables were used in an attempt to provide nicer
%alignment of columns, but this is better
%achieved by importing the tab delimited text file into Excel
%defaultSpacer=20;
%useSpacer = defaultSpacer;


%set the report fieldnames to the first cell (which is
%itself a cell array of strings. If only printing one table
%then allow bypassing the cell nesting.
if(iscell(params.report.tables.columns{1}))
  displayColumns = params.report.tables.columns{1};
else
  displayColumns = params.report.tables.columns;
end

%the only allowed argument that is not a fieldname is
%'answer_format', which lists the answer type. If used, this argument
%will first check enum_values, so we'll treat this as an
%'enum_values' argument
[answerFormatIsListed,answerFormatLoc] = ismember('answer_format',displayColumns);
if(answerFormatIsListed)
  displayColumns{answerFormatLoc} = 'enum_values';
end

numColumns = length(displayColumns);
numRecords = length(data_st);

%The first column of this 2-D cell array stores actual fieldnames
%(also used as tags to the 'columns' report field). The second
%column stores the string that is used to reformat the column name
%for nicer printout.
headingReformat = { 'experiment_id','Experiment ID';
		    'start_date','Start Date';
		    'experiment_title','Experiment Title';
		    'experiment_description','Experiment Description';
		    'response_table','Response Table';
		    'irb_id','IRB ID';
		    'end_date','End Date';
		    'language','Language';
		    'locked','Locked';
		    'form_id','Form ID';
		    'form_name','Form Name';
		    'form_category','Form Category';
		    'header','Header';
		    'footer','Footer';
		    'version','Version';
		    'locked','Locked';
		    'condition','Condition';
		    'condition_matlab','Condition Matlab';
		    'visit_once','Visit Once';
		    'compqid','Composite Question ID';
		    'question_id','Question ID';
		    'question_text','Question Text';
		    'question_category','Question Category';
		    'heading_format','Heading Format';
		    'locked','Locked';
		    'data_format_id','Data Format ID';
		    'type','Type';
		    'enum_values','Answer Format';
		    'subquestion','Subquestion';
		    'heading','Heading';
		    'range','Range';
		    'default','Default';
		    'html_field_type','HTML Field Type';
		    'required','Required'};
		    

[tf,colRftIdx] = ismember(displayColumns,{headingReformat{:,1}}); 

for iColumn = 1:numColumns
  columnName = displayColumns{iColumn};
  columnHeading = headingReformat{colRftIdx(iColumn),2};
  fprintf(fid,'%s',columnHeading);
  
  if(iColumn ~= numColumns)
    fprintf(fid,'\t');
  end
end

fprintf(fid,'\n');

for iRecord = 1:numRecords
  for iColumn = 1:numColumns
    
    %if the display column is not a valid field, then throw an
    %error message
    if(~isfield(data_st(iRecord),displayColumns{iColumn}))
      error(sprintf('Field %s is invalid for specified table structure', displayColumns{iColumn}));
    end
    
    %get the field value. This will need to be reformatted during
    %printout, depending on its type
    fieldVal = getfield(data_st(iRecord),displayColumns{iColumn});
    %get rid of any carriage return, line feed, and tab characters
    %in the field for clean tab delimited output
    if (ischar(fieldVal) | iscell(fieldVal))
      fieldVal = strrep(fieldVal,'\r','');
      fieldVal = strrep(fieldVal,'\t','');
      fieldVal = strrep(fieldVal,'\n','');
    end
    
    %check the field value type and format accordingly
    if(isnumeric(fieldVal))
      %spacer = useSpacer-length(headingReformat{colRftIdx(iColumn),2})+1;
      spacer = [];
      fmtString = ['%' num2str(spacer) 'd'];
      fprintf(fid,fmtString,fieldVal);
      
      
    elseif(strcmp(displayColumns{iColumn},'enum_values'))
      
      
      if(~isempty(fieldVal))
	for iCell = 1:length(fieldVal)
	  if(iCell == 1)
	    %spacer = useSpacer-length(headingReformat{colRftIdx(iColumn),2})+1;
	    spacer = [];
	    fmtString = ['%' num2str(spacer) 'd=%s,'];
	  elseif(iCell == length(fieldVal))
	    fmtString = '%d=%s';
	  else
	    fmtString = '%d=%s, ';
	  end
	  
	  fprintf(fid,fmtString,iCell,fieldVal{iCell});
	end
      
      else
	
	dataType = getfield(data_st(iRecord),'type');
	switch dataType
	 case {'int16','int32','int64'}
	  fieldReport = 'Any Integer';
	 case 'double'
	  fieldReport = 'Any Double';
	 case {'varchar','text'}
	  fieldReport = 'Any String';
	 case 'date'
	  fieldReport = 'Any valid date';
	 case 'year'
	  fieldReport = 'Any valid year';
	 case 'null'
	  fieldReport = 'None';
	end
	
	%spacer = useSpacer-length(headingReformat{colRftIdx(iColumn),2})+1;
	spacer = [];
	fmtString = ['%' num2str(spacer) 's\t'];
 
	fprintf(fid,fmtString,fieldReport);
	
      end
    
    elseif(ischar(fieldVal))
      
      %remove CR, LF, and TABS from the string (could not find a
      %way to use string matching functions for this)
      fieldVal(find(fieldVal == CR_ASCII | fieldVal == LF_ASCII | fieldVal == TAB_ASCII)) = [];
      
      %spacer = useSpacer-length(headingReformat{colRftIdx(iColumn),2})+length(fieldVal);
      spacer = [];
      fmtString = ['%' num2str(spacer) 's'];
      fprintf(fid,fmtString,fieldVal);

    
    end
    
    if(iColumn ~= numColumns)
      fprintf(fid,'\t');
    end
    
  end
  
  fprintf(fid,'\n');
end




return
