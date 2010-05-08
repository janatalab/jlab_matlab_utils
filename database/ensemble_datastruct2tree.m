function outtree = ensemble_datastruct2tree(dataStruct,params)
% Converts from an Ensemble data structure to a tree structure
% in the same format that mysql_extract_metadata uses.
%

% 16 March, 2007 Stefan Tomic
% 28 April, 2007 Petr Janata - changed tree to outtree, and added check to
%                make sure outtree isn't empty prior to attempting the call
%                to convert_structarray
%  5 July, 2007  Stefan Tomic - allowing empty data fields and
%                               differing data lengths. No longer
%                               performing conversion to array of structs
% 10/19/07  PJ Fixed trapping of empty values
% 2/22/10   ST Fixed conversion of regular (non-cell) struct arrays
  
outtree = [];

fnames = dataStruct.vars;

for ifld = 1:length(fnames)
  
  if(ifld > length(dataStruct.data))
    dataVals = [];
  else
    dataVals = dataStruct.data{ifld};
  end
  
  %dealing with whether a child structure is embedded in a cell or not
  if(iscell(dataVals) & (length(dataVals) >= 1) & isstruct(dataVals{1}))
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
outtree.name = dataStruct.name;
outtree.type = dataStruct.type;
outtree.report = dataStruct.report;
outtree.meta = dataStruct.meta;