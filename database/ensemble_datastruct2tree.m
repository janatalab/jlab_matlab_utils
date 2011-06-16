function outtree = ensemble_datastruct2tree(dataStruct,params)
% Converts from an Ensemble data structure to a tree structure in the same format that mysql_extract_metadata uses.
%
% outtree = ensemble_datastruct2tree(dataStruct,params);
%
% WARNING: If one of the data fields in the data structure has the same name
% as one of the datastruct fields, e.g. name, type, meta, the data will be
% overwritten unless a method of conflict resolution other than
% 'datastruct_rules' is specified in params.conflict.(fieldname), where
% fieldname is the name of the field in question.
%
% The reason that the datastruct field names are tacked on to the output
% structure is to preserve back-and-forth convertability between the tree
% and datastruct format. This is not optimal.
%
% If you know you don't need to preserve convertability, then go ahead and
% set a non-empty conflict resolution rule

% 16 March, 2007 Stefan Tomic
% 28 April, 2007 Petr Janata - changed tree to outtree, and added check to
%                make sure outtree isn't empty prior to attempting the call
%                to convert_structarray
%  5 July, 2007  Stefan Tomic - allowing empty data fields and
%                               differing data lengths. No longer
%                               performing conversion to array of structs
% 10/19/07  PJ Fixed trapping of empty values
% 2/22/10   ST Fixed conversion of regular (non-cell) struct arrays
% 16Jun2011 PJ Added conflict resolution
  
outtree = [];

defaultResolution = 'datastruct_rules';

fnames = dataStruct.vars;

for ifld = 1:length(fnames)
	
	if(ifld > length(dataStruct.data))
		dataVals = [];
	else
		dataVals = dataStruct.data{ifld};
	end
	
	%dealing with whether a child structure is embedded in a cell or not
	if(iscell(dataVals) && (length(dataVals) >= 1) && isstruct(dataVals{1}))
		for idx = 1:length(dataVals)
			if(isfield(dataVals{idx},'vars'))
				%perform ensemble_datastruct2tree if this is an ensemble
				%data struct
				dataVals{idx} = ensemble_datastruct2tree(dataVals{idx});
			end
		end
	end
	
	%deal with struct not in a cell
	if(isstruct(dataVals))
		clear convertedData;
		for idx = 1:length(dataVals)
			if(isfield(dataVals(idx),'vars')) && ~isempty(dataVals(idx).vars)
				%perform ensemble_datastruct2tree if this is an ensemble
				%data struct
				convertedData(idx) = ensemble_datastruct2tree(dataVals(idx));
			else
				convertedData(idx) = dataVals(idx);
			end
			
		end
		clear dataVals
		dataVals = convertedData;
	end
	
  outtree.(fnames{ifld}) = dataVals;
  
end

%copy ensemble data struct fields (name,type,report,meta)

% Check for naming conflicts, i.e. if a fieldname matches one of the
% ensemble struct data field names
ensemble_struct_names = {'name','type','report','meta'};
conflict_mask = ismember(ensemble_struct_names, fnames);
conflicts = ensemble_struct_names(conflict_mask);

if ~isempty(conflicts)
	warning('Found %d conflicting field names: %s\n', length(conflicts), cell2str(conflicts,','));
	
	% Try conflict resolution
	for ic = 1:length(conflicts)
		currConflict = conflicts{ic};
		
		% See if we've specified a resolutionMethod
		if exist('params','var') && isfield(params, 'conflict') && isfield(params.conflict, currConflict)
			resolutionMethod = params.conflict.(currConflict);
		else
			resolutionMethod = '';
		end
		
		% Report on resolution method and indicate use of default resolution
		% otherwise
		if ~isempty(resolutionMethod)
			fprintf('Resolving conflict for "%s" using method "%s"\n', currConflict, resolutionMethod);
		else
			fprintf('No resolution method specified for conflict "%s". using "%s"\n', currConflict, defaultResolution);
			resolutionMethod = defaultResolution;
		end
		
		% Implement resolution method
		switch resolutionMethod
			case 'datastruct_rules'
				outtree.(currConflict) = currConflict;
			otherwise
				% Leave the field as is
		end
	end
end

% Copy all non-conflict fields
for fld = ensemble_struct_names(~conflict_mask)
	outtree.(fld{1}) = fld{1};
end

return