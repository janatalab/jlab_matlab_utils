function s_out = convert_structarray(s_in)
%
% Accepts a struct with fields containing vectors of values, where
% the order of the elements of each vector link the values in the 
% fields to one another. This function converts this type of struct
% to an array of structs, where each instance in the array contains 
% a single value from the vectors contained in the original fields.
% 
% i.e. - a struct such as the following:
%
% myStruct = 
%       field1: [100 x 1 double]
%       field2: [100 x 1 double]
%       field3: {100 x 1 cell}
%
%
% is converted to the following:
%
% myConvertedStruct = 
%       1x100 struct array with fields:
%             field1
%             field2
%             field3
%
% 16 Jan, 2007 - First Version, S.T.



if(~isstruct(s_in) | (length(s_in) > 1))
  error('input argument must be a struct of length 1');
end


fldnames = fieldnames(s_in);

%first check the the vector lengths of all fields are the same
%before proceeding
for ifld = 1:length(fldnames)
  numElements(ifld) = length(getfield(s_in,fldnames{ifld}));
end

if(sum(diff(numElements)) ~= 0)
  error(['The lengths of all fields in the struct must be the same' ...
	 ' to use the function']);
end

newStruct = mkstruct(fldnames);

s_out = repmat(newStruct,1,numElements(1));


for ifld = 1:length(fldnames)
  valVector = getfield(s_in,fldnames{ifld});
  for iElement = 1:numElements
    if(iscell(valVector(iElement)))
      val = valVector{iElement};
    else
      val = valVector(iElement);
    end
    s_out(iElement) = setfield(s_out(iElement),fldnames{ifld},val);
  end
end
