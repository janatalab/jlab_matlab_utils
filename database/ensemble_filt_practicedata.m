function expinfo_practiceFilt = ensemble_filt_practicedata(expinfo,params)

% Removes practice trial data from stim_meta, response_data, and stim_Xattribute 
% structs produced by ensemble_load_expinfo.m. Output is expinfo struct
% with practice-related data filtered out. This is useful when some of the 
% trials presented in an Ensemble-controlled experiment are practice trials
% and one does not wish to include those trials in subsequent analyses.
%
% INPUT
%   expinfo (as produced by ensemble_load_expinfo.m)
%   params.practice_attrib -- attribute name associated with practice trials in Ensemble
%
% Oct 19, 2012 BH

[practice_id,] = mysql_get_stim_by_attribute(...
  'attrib_name', params.practice_attrib, ...
  'mysql', params.mysql);

expinfoCols = set_var_col_const(expinfo.vars);
expinfo_practiceFilt = expinfo;

% remove practice trial data from stimulus_metadata
stim_meta_st = expinfo.data{expinfoCols.stimulus_metadata};
stim_meta_cols = set_var_col_const(stim_meta_st.vars);
stim_ids = stim_meta_st.data{stim_meta_cols.stimulus_id};
stim_id_practice_mask = ismember(stim_ids, practice_id{1});
stim_meta_practiceFilt = stim_meta_st;
for icell = 1:length(stim_meta_st.data)
    stim_meta_practiceFilt.data{icell}(stim_id_practice_mask) = [];
end
expinfo_practiceFilt.data{expinfoCols.stimulus_metadata} = stim_meta_practiceFilt;

% remove practice trial data from response_data
resp_st = expinfo.data{expinfoCols.response_data};
resp_st_cols = set_var_col_const(resp_st.vars);
resp_practice_mask = ismember(resp_st.data{resp_st_cols.stimulus_id},practice_id{1});
resp_data_practiceFilt = resp_st;
for iRespCell = 1:length(resp_st.data)
    resp_data_practiceFilt.data{iRespCell}(resp_practice_mask) = [];
end
expinfo_practiceFilt.data{expinfoCols.response_data} = resp_data_practiceFilt;

% remove practice trial data from stimulus_x_attribute
stimXatt_st = expinfo.data{expinfoCols.stimulus_x_attribute_metadata};
stimXatt_cols = set_var_col_const(stimXatt_st.vars);
stimXatt_practice_mask = ismember(stimXatt_st.data{stimXatt_cols.stimulus_id},practice_id{1});
stimXatt_practiceFilt = stimXatt_st;
for iStimAttCell = 1:length(stimXatt_st.data)
    stimXatt_practiceFilt.data{iStimAttCell}(stimXatt_practice_mask) = [];
end
expinfo_practiceFilt.data{expinfoCols.stimulus_x_attribute_metadata} = stimXatt_practiceFilt;

return