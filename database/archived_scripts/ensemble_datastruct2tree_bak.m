function outtree = ensemble_datastruct2tree(dataStruct,params)
%
%
% converts from an Ensemble data structure to a tree structure in
% the same format that mysql_extract_metadata uses.
%

% 16 March, 2007 Stefan Tomic
% 28 April, 2007 Petr Janata - changed tree to outtree, and added check to
%                make sure outtree isn't empty prior to attempting the call
%                to convert_structarray
%              

outtree = [];

fnames = dataStruct.vars;

for ifld = 1:length(fnames)
  
  dataVals = dataStruct.data{ifld};
  
  if(iscell(dataVals) & isstruct(dataVals{1}))
    
    
    for idx = 1:length(dataVals)

      dataVals{idx} = ensemble_datastruct2tree(dataVals{idx});      
      
    end
    
    
  end
  
  outtree.(fnames{ifld}) = dataVals;
  
end

if ~isempty(outtree)
    outtree = convert_structarray(outtree);
end


  
