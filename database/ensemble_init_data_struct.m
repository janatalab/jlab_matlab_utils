function data_st = ensemble_init_data_struct
% Initializes a data structure used by many of the Ensemble analysis functions.
%
% data_st = ensemble_init_data_struct;
%
% Fields are:
%    .type - This is a string that identifies the data structure
%    .vars - A cell array of strings that identifies each of the elements of
%            the data cell array. This variable is used for indexing into the
%            data array and also by ensemble_filter().
%    .data - A cell array containing either data structures or data
%    .report - A structure that contains information pertaining to reporting
%            tables and figures related to this type of data
%            structure. Analysis functions that have reporting sub-functions
%            should register function handles within this structure.
%    .meta - Structure for holding metadata pertaining to this type of data
%            structure.  This might include the parameter structure that was
%            passed to the function that generated this structure, or a
%            question information structure, etc.
%
%
% Copyright (c) 2007-2012 The Regents of the University of California
% All Rights Reserved
%
% 02/03/07 Petr Janata

data_st.name = '';
data_st.type = '';
data_st.vars = {};
data_st.data = {};
data_st.report = struct();
data_st.meta = struct();;

