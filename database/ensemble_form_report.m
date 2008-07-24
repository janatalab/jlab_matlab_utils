function expDataStruct = ensemble_form_report(expDataStruct,params)
% Prints out a report of forms and questions used in an experiment
% This is sort of a wrapper for ensemble_print_metadata, in that it 
% filters out the forms you want and restructures the data in the
% tree structure that ensemble_print_metadata can read
%
%

fnames = expDataStruct.vars;
if (~ismember('form',fnames))
  error('This data structure does not contain form data');
end

metaFieldConst = set_var_col_const(fnames);
formDataStruct = expDataStruct.data{metaFieldConst.form};

formFilter = params.filt;


for(fIdx = 1:length(formDataStruct))
  formDataStruct{fIdx} = ensemble_filter(formDataStruct{fIdx},formFilter);
end

%assign filtered forms back to expDataStruct
expDataStruct.data{metaFieldConst.form} = formDataStruct;


expTree = ensemble_datastruct2tree(expDataStruct);


printParams.report.tables.write2file = 1;
printParams.report.tables.fname      = '/tmp/mymeta.csv';
printParams.report.tables.columns = {{},...
		    {'form_name'},...
		    {'question_text','enum_values'}};
ensemble_print_metadata(expTree,printParams);
