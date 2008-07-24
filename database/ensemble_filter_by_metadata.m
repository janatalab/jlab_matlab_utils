function data_st = ensemble_filter_by_metadata(data_st,params)
% Filters a given data set by given metadata.
% accepts a data struct and corresponding metadata, which are
% stored in data_st as a cell array. Accepts filt params as
% ensemble_filter
%
% 14, Feb. 2007 S.T.

filt = params.filt;

inc_exc = fieldnames(filt);


all_any = fieldnames(filt.(inc_exc{1}));

crit_fld = fieldnames(filt.(inc_exc{1}).(all_any{1}));

%go through each data structure and filter based on crit field (if
%found in the data structure
for ist = 1:length(data_st)
  idxCritFld = strmatch(crit_fld{1},data_st{ist}.vars,'exact');
  if(~isempty(idxCritFld))
    data_st = ensemble_filter(data_st,params);
  end
end


%find the fields that the data structures have in common and take
%the intersection of the records with common values
%NOTE: RIGHT NOW THIS IS ONLY DEALING WITH A PAIR OF DATA
%STRUCTURES. IT WILL NEED TO BE EXPANDED TO DEAL WITH MORE THAN 2.
[tf,loc] = ismember(data_st{1}.vars,data_st{2}.vars);

for colIdx = 1:length(loc)
  
  if(loc(colIdx) ~= 0)
    
    fldIdx(1) = colIdx
    fldIdx(2) = loc(colIdx);
    
    crit_vals = intersect(data_st{1}.data{fldIdx(1)},data_st{2}.data{fldIdx(2)});
    %if the column has more values than crit_vals, then filter this
    %data_st
    
    for dataStIdx = 1:2
      if(~isempty(setdiff(data_st{dataStIdx}.data{fldIdx(dataStIdx)},crit_vals)))
      
	critFiltParams.filt.include.all.(crit_fld{1}) = crit_vals;
	
	data_st{dataStIdx} = ensemble_filter(data_st{dataStIdx},critFiltParams);
      
      end
    end
  end
end
