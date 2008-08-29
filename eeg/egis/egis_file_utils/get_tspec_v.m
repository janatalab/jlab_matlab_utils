function spec_value = get_tspec_v(cellhdrs,cell,trial, feild)
%spec_value=get_tspec(cell_hdrs,cell,trial, feild)
%
%Returns a value which is an EGIS trial-specific from a matrix cell_hdrs.
%This matrix is that returned as the second value of the function
%rd_ses_hdr_v.
%
%See Also:rd_ses_hdr_v, get_tspecs_v.

tspec_offset = 5 + ((trial-1)*(cellhdr(cell,5)/2)) + feild; 
spec_value = cell_hdrs(cell,tspec_offset);
