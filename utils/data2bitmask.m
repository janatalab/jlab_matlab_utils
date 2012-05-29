function bitmask = data2bitmask(data, max_bit)
% data2bitmask - converts a matrix of numeric data values into bitmasks that
% correspond to those values.
%
% maxbit - the maximum bit than can be set.  This is important insofar that it
% pads the mask with extra columns in the event that none of the data values
% set the highest bit
%
% bitmask is a logical vector

% 09/09/05 Petr Janata

% convert the number to a binary representation as a string of zeros and ones
bin = dec2bin(data);

% convert the ascii characters to their numeric values and subtract the numeric
% value for '0'
bitmask  = double(bin)-double('0');

% since the most significant bit is currently on the left, we want to flip the
% matrix left to right so that the column number corresponds to the bit number
bitmask = fliplr(bitmask);

% if a maximum bit has been specified, make sure that the number of columns in
% our mask corresponds to the desired number of bits
if exist('max_bit','var') && max_bit > size(bitmask,2)
  bitmask(end,max_bit) = 0;
end

bitmask = logical(bitmask);

return
