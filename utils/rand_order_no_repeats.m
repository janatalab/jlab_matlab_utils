function [output_vector, vector_ind] = rand_order_no_repeats(input_vector)
% randomizes the order of the input_vector and attempts to not allow any one
% value to be repeated in consecutive order (assuming that there are repeats
% in the values). Note that this won't be possible if the number of occurrences
% of any one value is greater than round(length(vector_ind)/2)

rand('state',sum(100*clock));

rand_input_vector = input_vector(randperm(length(input_vector)));
unique_values = unique(rand_input_vector);

for(idx = 1:length(unique_values))
  vals(idx).value = unique_values(idx);
  vals(idx).order = find(rand_input_vector == vals.value);
  vals(idx).repeats = find(vals(idx).order == vals(idx).order+1);
end



%now do the work of checking for consecutive values and reorder appropriately
for(unique_idx = 1:length(unique_values))
  
  if(~isempty(vals(unique_idx).repeats))
    
    for(swap_idx = setdiff([1:length(unique_values)],unique_idx))
      
      for(repeat_idx = vals(unique_idx).repeats)
	
	if(
	
	
	
	
      end
      
      
      
      
      
    end
    
    
    
  
  end
  
end