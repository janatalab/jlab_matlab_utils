function outData = ensemble_export_respnstim(inData,params)

% tabulates subject-level and stimulus-based response-level data, exports one row per stim
%
% function outData = ensemble_export_respnstim(inData,params)
%
% This function takes given response table data, extracts data that are
% subject-level data (such as scores on a personality scale), extracts
% stimulus-level data (such as the responses made to each stimulus) and
% returns data in a form where there is one row per stimulus presented per
% subject, and any subject-level scores are applied to each row for that
% subject. Data can be saved to a demilited file (comma default, other
% delimiters can be specified), and an optional SAS script can be generated
% to import the data into a SAS dataset.
% 
% INPUT/PARAMS
% 
%  inData should be a cell array of structs, containing the following
%    possible data:
%    - the 'response_data' struct from ensemble_load_expinfo - REQUIRED
%    - the 'stimulus_metadata' struct from ensemble_load_expinfo - OPTIONAL
%    - the 'subject_info' struct from ensemble_load_expinfo - OPTIONAL
%    - the results of any function that has calculated subject-level
%    constants and contains a subject_id column. currently, all vars other
%    than subject_id that are passed on will be included in the output data
%    as constants for that subject.
% 
%  params.filt - any filtering criteria found here are applied at once to
%    the response_data struct. This can be utilized to include or exclude
%    forms, questions, subquestions, sessions, subjects, etc (AND compqids)
%    NOTE: make sure not to limit your filtering parameters so that only
%    questions for a given stimulus or form are included, IF you are
%    passing that data onto a function that is called within this script
%    (through the params.export.by_subject.(fieldname) or
%    .by_stimulus.(fieldname) mechanisms) for subject or stimulus-level
%    constants
% 
%  params.export.by_subject.(fieldname) - fieldname can either be the
%    analysis struct name of a dataset within inData, or the name of a
%    function to be executed to return a dataset, containing a subject_id
%    variable and data that you would like to use as subject-level constant
%    data. - OPTIONAL
%    NOTE: if any given dataset contains more than one row per
%    subject_id, a warning will be issued, and this script will output only
%    the first observation for that subject_id
%  params.export.by_subject.(fieldname).vars - a cell array of strings of
%    vars to export, from the dataset provided by (fieldname) - OPTIONAL
%    NOTE: if vars in this cell array do not exist within the output of
%    (fieldname), then all vars in (fieldname) will be exported as
%    subject-level constants
% 
%  params.export.by_stimulus.(fieldname) - similar to
%    params.export.by_subject.(fieldname), except for stimulus-level
%    metadata ... the return of this function handle or inData struct will
%    be treated as stimulus-level constants - OPTIONAL
%  params.export.by_stimulus.(fieldname).vars - similar to
%    params.export.by_subject.(fieldname).vars, except for the output of
%    params.export.by_stimulus.(fieldname) - OPTIONAL
% 
%  params.export.non_num_delim - if present, the value of this param is
%    placed before and after non-numeric cell values when exporting (i.e.
%    single or double quotes around non-numeric cells)
%  params.export.delimiter - character to delimit cell values. if none is
%    provided, a tab is used
%  params.export.var_name_map - a struct whose fieldnames are the
%    calculated fieldnames from question/subquestion ids. These fieldnames
%    will be located within the outData.vars array and converted to their
%    values within params.export.var_name_map. For instance, if
%    params.export.var_name_map.s123_01 = 'renamedq', and if question 123,
%    subquestion 1 exists within the dataset, the varname for this question
%    in the output will be renamed to 'renamedq'. Original var names are
%    returned in the outData.origvars field.
% 
%  params for ensemble_init_fid
%    params.write2file, params.print, params.fname, params.filemode
%  params.sas.fname - REQUIRED if you want to automatically generate a SAS
%    import file for your data.
%  params.sas.options - an optional 'options' string for the beginning of
%    your sas file
%  params.sas.libsave - if you set this = 1, and you have defined
%    params.sas.libpath and params.sas.libname, it will include the command
%    to save your imported data to a SAS data library.
%  params.sas.libpath - the path to your SAS data library directory
%  params.sas.libname - the filename that will be given to your new SAS
%    data library, if you have chosen to save this dataset as a library
%
% OUTPUT
%  outData - contains one row per stimulus per subject, one column for
%    subject_id, stimulus_id, and trial_id, one column per subject-level
%    score calculated, one column per stimulus_metadata var that was
%    included in the inData struct, and one column per stimulus-level qnum.
%  outData.origvars - the original variable names, most of which will be
%    constructed from question and subquestion numbers. This is provided as
%    a reference for re-named variables (see params.
% 
% NOTE: the function gets it's list of valid subject_ids from
% response_data, so if this is not required, the function will return
% empty. But ... if you provide response_data, and then do not specify any
% params.export.by_stimulus fields, then it will only return subject-level
% data. this is a bit messy, a better way to do it would be to glean
% subject ids from the 'subject_info' struct, and make it a hard
% requirement, but that is not implemented yet ... it was coded this way at
% first since this is supposed to be mainly a function with a
% stimulus-response level focus, but there is room for improvement here so
% I will say that this warrants a FIXME
% 
% NOTE: if a given subject-level score variable contains more than one of
% any given score for a subject, the first of each will be used. also, if
% more than one iteration of a qnum exists for any given question/stimulus
% pair, the first of the list will be reported for that stimulus
% 
% NOTE: this script assumes that a given question will only be asked once
% for a given stimulus_id within a given subject. If the same compqid is
% found more than once for a given stimulus within a given subject's data,
% the first response will be returned and the subsequent responses will be
% ignored.
% 
% FIXME - need to extract/integrate data from the misc_info column of a
% given response table, include in the row for each stimulus ... up to this
% point, I have dealt with data stored in misc info by writing a separate
% function to extract that data and add it as another variable in the
% response_data struct. This might be a good way to deal with misc_info dat
% moving forward.
%
% FB 12/17/07 - started scripting
% FB 10/25/09 - added support for multiple trials with the same stimulus.
% If a stimulus is presented multiple times, the trial number should be
% incremented, and if this is done, export_respnstim will use the trial
% number to dissociate responses to the same stimulus at different
% presentations.
% FB 01/22/10 - added 'response_order' to the output for 'by_stimulus'
% output ... now outputs the response order of the first response for a
% given session/subject/stimulus/trial. In this way, if you sort your data
% by stimulus, or any other variable, you can later retrieve the order of
% presentation of stimuli by sorting by response_order.

% % initialize output data struct
outData = ensemble_init_data_struct;
outData.type = 'resp_x_stim'; % unique analysis type
outData.vars = {}; % names of columns that will end up inas outData.data
outData.datatype = {}; % tracks non-numeric data, for SAS output
if isfield(params,'outDataName') && ~isempty(params.outDataName)
    outData.name = params.outDataName;
else
    outData.name = 'resp_x_stim';
end;

% % check for export criteria
if ~isfield(params,'export')
    % if no export criteria are provided, we can not process export data
    warninig('no export params specified');
    return
end

% % check inData, convert to cell if needed
if ~iscell(inData)
    inData = {inData};
end

% % get response data, filter
rdSearchCrit.name = 'response_data';
respIdx = ensemble_find_analysis_struct(inData,rdSearchCrit);
if respIdx
    rsData = inData{respIdx};
    % check for, generate compqids
    rsData = ensemble_check_compqid(rsData);
else
    % if no response data are present, you can't export them
    warning('no response_data struct was found');
    return
end

% apply filtering criteria, if present
if isfield(params,'filt')
  rsData = ensemble_filter(rsData,params.filt);
end

rsCols = set_var_col_const(rsData.vars);

% % initialize subject-level data structure
% sNum is an auto-incrementing subject number within this dataset
sub_st.vars={'subject_id','sNum'};
sub_st.data{1} = unique(rsData.data{rsCols.subject_id});
subs = sub_st.data{1};
nsub = length(sub_st.data{1});
sub_st.data{2} = 1:nsub;
% sub_st.datatype tracks data type
%   (s)tring
%   (n)umeric
%   (l)ogical
sub_st.datatype = {'s','n'};

% % extract/calculate subject-level scores
if isfield(params.export,'by_subject')

    % % cycle through params.export.by_subject fieldnames
    % % if struct is found in inData, load that data (expects fldname)
    % % if struct is NOT found in inData, try to excute it as a function
    % % if no subject_id var in output or struct, warn & continue
    % % if include.(fldname).vars, include only those vars
    subinc = params.export.by_subject;
    subfldnames = fieldnames(subinc);
    for ifld=1:length(subfldnames)
        lsData = {};
        lfld = subfldnames(ifld);
        if iscell(lfld), lfld = lfld{1}; end
        lSchCrit.name = lfld;
        lIdx = ensemble_find_analysis_struct(inData,lSchCrit);
        if lIdx
            lsData = inData{lIdx};
        else
            lfh = str2func(lfld);
            try lsData = lfh(rsData,params);
            catch
                lsData = {};
            end;
        end 
        if ~isempty(lsData) && isstruct(lsData)
            lsCols = set_var_col_const(lsData.vars);

            % sanity check, make sure one row per subject
            if isfield(lsCols,'subject_id')
                usubs = length(unique(lsData.data{lsCols.subject_id}));
                if usubs
                    nrows = length(lsData.data{lsCols.subject_id});
                    if nrows ~= usubs
                        warning(['subject-level data (%s) did not have '...
                            'one row per subject'],lfld)
                    end
                end
            else
                warning(['subject-level data (%s) did not have '...
                    'a val named ''subject_id'''],lfld)
                continue
            end
            % if 'vars' field found under params.export.by_subject.(lfld)
            % then mask out all other vars & data
            if isfield(subinc.(lfld),'vars') && iscell(subinc.(lfld).vars)
                xvars = subinc.(lfld).vars;
                if ~ismember(lsData.vars,xvars)
                    warning(['vars defined for %s can not be found, '...
                        'exporting all vars for this subject-level '...
                        'data source',lfld]);
                else
                    if ~ismember(xvars,'subject_id')
                        xvars = [xvars 'subject_id'];
                    end
                    xmask = ismember(lsData.vars,xvars);
                    lsData.vars = lsData.vars(xmask);
                    lsData.data = lsData.data(xmask);
                end
            end
            lsMask = ~ismember(lsData.vars,{'subject_id','session_id'});
            sub_st.vars = [sub_st.vars lsData.vars{lsMask}];
            lsIdxs = find(lsMask);
            % check, track, data type (s)tring or (n)umeric
            for iidx = 1:length(lsIdxs)
                 if isnumeric(lsData.data{lsIdxs(iidx)})
                     dtype = 'n';
                 elseif islogical(lsData.data{lsIdxs(iidx)})
                     dtype = 'l';
                 else
                     dtype = 's';
                 end
                 sub_st.datatype = [sub_st.datatype dtype];
            end
            vidxs = find(ismember(sub_st.vars,{lsData.vars{lsMask}}));
            nvars = length(vidxs);
            for ivar = 1:nvars
                vi = vidxs(ivar);
                if sub_st.datatype{vi} == 's'
                    sub_st.data{vi} = cell(nsub,1);
                else
                    sub_st.data{vi} = zeros(nsub,1);
                end
            end

            % parse rows
            for isub = 1:nsub
                sidx = find(ismember(lsData.data{lsCols.subject_id},subs{isub}));
                if length(sidx) > 1
                  sidx = sidx(1);
                elseif isempty(sidx)
                  continue
                end
                for ivar = 1:nvars
                    vi = vidxs(ivar);
                    li = find(ismember(lsData.vars,sub_st.vars{vi}));
                    if sub_st.datatype{vi} == 's'
                        svdata = lsData.data{li}{sidx};
                    else
                        svdata = lsData.data{li}(sidx);
                    end
                    svdata = sanitize_cell_value(svdata,sub_st.datatype{vi});
                    if sub_st.datatype{vi} == 's'
                        sub_st.data{vi}{isub} = svdata;
                    else
                        sub_st.data{vi}(isub) = svdata;
                    end
                end
            end
        else
            warning('couldn''nt find inData struct or function handle for %s',...
                lfld);
        end % for ~isempty(lsData
    end % for ifld=1:length(subfldnames
end % if isfield(params.export,'by_subject'

% initalize vars
outData.vars = sub_st.vars;
outData.datatype = sub_st.datatype;

% % extract/calculate stim-level scores
% filter by each stim-level qnum, loop over each subject, data2enum and
% store scores
if isfield(params.export,'by_stimulus')

    % add var for trial_id, stim repetition number, and response order
    outData.vars = [outData.vars 'trial_id' 'stim_rep' 'response_order'];
    outData.datatype = [outData.datatype 'n' 'n' 'n'];
    
    % % initialize stimulus-level data structure
    smData.vars={'stimulus_id'};
    smData.datatype = {'n'}; % tracks data type, (s)tr, (n)um, (l)ogical
    smData.data = {[]};
    smCols = set_var_col_const(smData.vars);
    smIdxs = [];
    stims = [];

    if isstruct(params.export.by_stimulus)
        stiminc = params.export.by_stimulus;
        stimfldnames = fieldnames(stiminc);
        for ifld=1:length(stimfldnames)
            sfld = stimfldnames(ifld);
            if iscell(sfld)
                sfld = sfld{1};
            end
            sSchCrit.name = sfld;
            sIdx = ensemble_find_analysis_struct(inData,sSchCrit);
            if sIdx
                lsData = inData{sIdx};
            else
                sfh = str2func(sfld);
                try lsData = sfh(rsData,params);
                catch
                    lsData = {};
                end;
            end

            if ~isempty(lsData) && isstruct(lsData)
                lsCols = set_var_col_const(lsData.vars);

                % sanity check, make sure one row per stimulus
                if isfield(lsCols,'stimulus_id')
                    ustim = unique(lsData.data{lsCols.stimulus_id});
                    if ~isnumeric(ustim)
                        warning(['stimulus_id colum for %s is non-numeric, '...
                            'please use only numeric stimulus_ids. %s will '...
                            'not be included as stim-level constant data'],...
                            sfld,sfld);
                        continue;
                    end
                    if ustim
                        nrows = length(lsData.data{lsCols.stimulus_id});
                        if nrows ~= length(ustim)
                            warning(['stimulus-level data (%s) did not have '...
                                'one row per stimulus'],sfld)
                        end
                    end
                    % add stimulus ids that aren't already in smData
                    newstim = setdiff(ustim,stims);
                    if newstim
                        smData.data{smCols.stimulus_id} = [stims newstim];
                        stims = smData.data{smCols.stimulus_id};
                        nstim = length(stims);
                        % initialize smData.data cols for each new row
                        for ivar = 2:length(smData.vars)
                            if smData.datatype{ivar} == 's'
                                smData.vars{ivar}{nstim} = [];
                            else
                                smData.data{ivar}(nstim) = 0;
                            end
                        end
                    end
                else
                    warning(['stimulus-level data (%s) did not have '...
                        'a val named ''stimulus_id'''],sfld)
                    continue
                end

                % if 'vars' field found under params.export.by_stimulus.(sfld)
                % then mask out all other vars & data
                if isfield(stiminc.(sfld),'vars') && iscell(stiminc.(sfld).vars)
                    xvars = stiminc.(sfld).vars;
                    if ~any(ismember(lsData.vars,xvars))
                        warning(['vars defined for %s can not be found, '...
                            'exporting all vars for this stimulus-level '...
                            'data source',sfld]);
                    else
                        if ~ismember(xvars,'stimulus_id')
                            xvars = [xvars 'stimulus_id'];
                        end
                        xmask = ismember(lsData.vars,xvars);
                        lsData.vars = lsData.vars(xmask);
                        lsData.data = lsData.data(xmask);
                    end
                end
                lsMask = ~ismember(lsData.vars,{'stimulus_id'});
                lsIdxs = find(lsMask);
                smData.vars = [smData.vars lsData.vars{lsMask}];
                smCols = set_var_col_const(smData.vars);
                smIdxs = find(~ismember(smData.vars,'stimulus_id'));
                % check, track, data type (s)tring or (n)umeric
                for iidx = 1:length(lsIdxs)
                     if isnumeric(lsData.data{lsIdxs(iidx)})
                         dtype = 'n';
                     elseif islogical(lsData.data{lsIdxs(iidx)})
                         dtype = 'l';
                     else
                         dtype = 's';
                     end
                     smData.datatype = [smData.datatype dtype];
                end
                vidxs = find(ismember(smData.vars,{lsData.vars{lsMask}}));
                nvars = length(vidxs);
                for ivar = 1:nvars
                    vi = vidxs(ivar);
                    if smData.datatype{vi} == 's'
                        smData.data{vi} = cell(nstim,1);
                    else
                        smData.data{vi} = zeros(nstim,1);
                    end
                end

                % parse rows
                for istim = 1:nstim
                    sidx = find(ismember(lsData.data{lsCols.stimulus_id},stims(istim)));
                    if length(sidx) > 1
                        sidx = sidx(1);
                    end
                    for ivar = 1:nvars
                        vi = vidxs(ivar);
                        li = find(ismember(lsData.vars,smData.vars{vi}));
                        if ~isempty(li)
                            if smData.datatype{vi} == 's'
                                svdata = lsData.data{li}{sidx};
                            else
                                svdata = lsData.data{li}(sidx);
                            end
                            svdata = sanitize_cell_value(svdata,smData.datatype{vi});
                            if smData.datatype{vi} == 's'
                                smData.data{vi}{istim} = svdata;
                            else
                                smData.data{vi}(istim) = svdata;
                            end % if smData.datatype{vi
                        end % if ~isempty(li
                    end % for ivar=1:
                end % for istim = 1:
            else 
                warning('couldn''nt find inData struct or function handle for %s',...
                    sfld);
            end % for ~isempty(lsData
        end % for ifld=1:length(stimfldnames
    end % isstruct(params.export.by_stimulus

    outData.vars = [outData.vars smData.vars];
    outData.datatype = [outData.datatype smData.datatype];

    cqids = {};
    ustims = {};
    % filter out non-stimulus response data
    ustims = unique(rsData.data{rsCols.stimulus_id});
    stimidxs = ~isnan(ustims);
    ustims = ustims(stimidxs);
    filt = {};
    filt.include.all.stimulus_id = ustims;
    aqData = ensemble_filter(rsData,filt);
    % get unique qnums - now compqids
    cqids.vars = {'compqid','question_id','subquestion',...
        'qinfo','scqid','bitmask'};
    cqCols = set_var_col_const(cqids.vars);
    cqScq = cqCols.scqid;
    cqCqi = cqCols.compqid;
    cqQid = cqCols.question_id;
    cqSub = cqCols.subquestion;
    cqQin = cqCols.qinfo;
    cqBit = cqCols.bitmask;
    ucqids = unique(aqData.data{rsCols.compqid});
    ucqididxs = find(~isnan(ucqids));
    nucqid = length(ucqididxs);
    for icqc = 1:length(cqids.vars)
        cqids.data{icqc} = cell(nucqid,1);
    end
    cqids.data{cqScq} = ucqids(ucqididxs);

    qnums = {};
    % parse out checkbox enums
    for icq = 1:length(cqids.data{cqScq})
        % check each cqid
        % could be converted to 'ensemble_compqid2dfid'
        [qidx] = find(ismember(rsData.data{rsCols.compqid},cqids.data{cqScq}(icq)));
        if length(qidx) > 1
            qidx = qidx(1);
        end
        cqids.data{cqQid}(icq) = num2cell(rsData.data{rsCols.question_id}(qidx));
        cqids.data{cqSub}(icq) = num2cell(rsData.data{rsCols.subquestion}(qidx));
        lcqid = cqids.data{cqScq}(icq);
        if iscell(lcqid)
            lcqid = lcqid{1};
        end
        cqids.data{cqCqi}(icq) = {sprintf('%1.2f',lcqid)};
        lqid = cqids.data{cqQid}(icq);
        if iscell(lqid)
            lqid = lqid{1};
        end
        lsqid = cqids.data{cqSub}(icq);
        if iscell(lsqid)
            lsqid = lsqid{1};
        end
        lqinfo = mysql_extract_metadata('table','question',...
            'question_id',lqid,'conn_id',params.mysql.conn_id);
        cqids.data{cqQin}(icq) = {lqinfo};
        uhft = {lqinfo.html_field_type};
        if strmatch('checkbox',uhft{lsqid});
            % get dfid for qid, expand qnums for this compqid
            ueval = {lqinfo.enum_values};
            for istr = 1:length(ueval{lsqid})
                lqnum = sprintf('%1.2f_c%02d',lcqid,istr);
                qnums = [qnums lqnum];
            end
        else
            qnums = [qnums cqids.data{cqCqi}(icq)];
        end
    end

    % add vars/init data cells
    qStructs = make_valid_struct_key(num2cell(qnums));
    numconst = length(outData.vars);
    for iadvar = numconst+1:numconst+length(qnums)
        outData.data{iadvar} = [];
        %%%% FIXME: assumes all responses will be numeric
        outData.datatype = [outData.datatype 'n'];
    end
    outData.vars = [outData.vars qStructs];
    outCols = set_var_col_const(outData.vars);    

    % a few constants
    scol = rsCols.stimulus_id;
    qcol = rsCols.question_id;
    ecol = rsCols.response_enum;
    tcol = rsCols.trial_id;
    rocol = rsCols.response_order;
    
    % sum(nstimpersub) = nrows
    nrows = 0;
    subData.vars = {'subid','locData','stim'};
    subCols = set_var_col_const(subData.vars);
    for isdv = 1:length(subData.vars)
        subData.data{isdv} = cell(nsub,1);
    end
    for isub = 1:nsub
        % get this subject's response data
        locFilt.include.all.subject_id = {subs{isub}};
        locData = ensemble_filter(rsData,locFilt);
        subData.data{subCols.subid}(isub) = {subs{isub}};
        subData.data{subCols.locData}(isub) = {locData};

        % get unique stim
        ustims = unique(locData.data{scol});
        stimidxs = ~isnan(ustims);
        stim = ustims(stimidxs);
        subData.data{subCols.stim}(isub) = {stim};

        nstim = length(stim);
        nrows = nrows + nstim;
    end

    for ivar = 1:length(outData.vars)
        if outData.datatype{ivar} == 's'
            outData.data{ivar} = cell(nrows,1);
        else
            outData.data{ivar} = zeros(nrows,1);
        end
    end
    
    % row counter
    outRow = 0;
    % % collate stim-level data
    for isub = 1:nsub
        subid   = subData.data{subCols.subid}(isub);
        locData = subData.data{subCols.locData}(isub);
        stim    = subData.data{subCols.stim}(isub);
        locData = locData{1};
        stim    = stim{1};
        nstim   = length(stim);

        for istim = 1:nstim
          stmFilt = {};
          stmFilt.include.all.stimulus_id = stim(istim);
          stmDataT = ensemble_filter(locData,stmFilt);
          uTrials = unique(stmDataT.data{tcol});
          nTrials = ~isnan(uTrials);
          nUt = sum(nTrials);
          if ~nUt, nUt = 1; end
          
          for in=1:nUt
            if any(nTrials)
              tFilt.include.all.trial_id = uTrials(in);
              stmData = ensemble_filter(stmDataT,tFilt);
            else
              stmData = stmDataT;
            end
            
            outRow = outRow + 1;
            
            % set id vars
            outData.data{outCols.subject_id}(outRow) = subid;
            outData.data{outCols.sNum}(outRow) = isub;
            outData.data{outCols.stimulus_id}(outRow) = stim(istim);
            outData.data{outCols.trial_id}(outRow) = uTrials(in);
            outData.data{outCols.stim_rep}(outRow) = in;
            
            % set subject-level scores
            % start at 3, since subject_id and sNum are always 1 and 2
            for icon = 3:length(sub_st.vars)
                outData.data{icon}(outRow) = sub_st.data{icon}(isub);
            end
            % set stimulus metadata
            if smIdxs
                stimidx = ismember(smData.data{smCols.stimulus_id},stim(istim));
                for ismd = 1:length(smIdxs)
                    smdidx = smIdxs(ismd);
                    smvar = smData.vars{smdidx};
                    lsmpt = smData.data{smdidx}(stimidx);
                    ltype = outData.datatype{outCols.(smvar)};
                    lsmpt = sanitize_cell_value(lsmpt,ltype);
                    if ltype == 's'
                        outData.data{outCols.(smvar)}{outRow} = lsmpt;
                    else
                        outData.data{outCols.(smvar)}(outRow) = lsmpt;
                    end
                end
            end

            % set response_order
            outData.data{outCols.response_order}(outRow) = stmData.data{rocol}(1);
            
            for iq = 1:length(qnums)
                qi = qnums(iq);
                if iscell(qi)
                    qi = qi{1};
                end
                % if not checkbox, go for it, if check box, parse out vals
                qcidx = regexp(qi,'_c\d{2}');
                if iscell(qcidx)
                    qcidx = qcidx{1};
                end
                qdata = {};
                if ~isempty(qcidx)
                    % data2bitmask, get items
                    lcqid = qi(1:qcidx-1);
                    lcqidx = find(ismember(cqids.data{cqCqi},lcqid));
                    lqinfo = cqids.data{cqQin}(lcqidx);
                    if (length(lqinfo)>1)
                        lqi = strmatch('checkbox',lqinfo.html_field_type);
                    else
                        lqi = 1;
                    end
                    evals = {lqinfo{lqi}.enum_values};
                    nenum = length(evals{lqi});
                    lqid  = {lqinfo{lqi}.question_id};
                    if iscell(lqid)
                        lqid = lqid{1};
                    end
                    bidx  = str2double(qi(qcidx+2:end));
                    lsqidx = find(ismember(stmData.data{qcol},lqid));
                    if length(lsqidx) > 1
                        lsubq = cqids.data{cqSub}(lcqidx);
                        lsubq = lsubq{1};
                        if length(lsqidx) < lsubq
                            lsqidx = lsqidx(1);
                        else
                            lsqidx = lsqidx(lsubq);
                        end
                    end
                    if lsqidx
                        qenum = stmData.data{ecol}(lsqidx);
                        if bidx <= nenum
                            if isempty(cqids.data{cqBit}{lcqidx})
                                qdata = data2bitmask(qenum,nenum);
                                cqids.data{cqBit}{lcqidx} = qdata;
                            else
                                qdata = cqids.data{cqBit}{lcqidx};
                            end
                            qdata = qdata(bidx);
                        else
                            qdata = '.';
                        end
                    else
                        qdata = '.';
                    end % if lsqid
                else
                    lqidx = find(ismember(cqids.data{cqCqi},qi));
                    lqid  = cqids.data{cqQid}(lqidx);
                    if iscell(lqid)
                        lqid = lqid{1};
                    end
                    [qidx] = find(ismember(stmData.data{qcol},lqid));
                    if length(qidx) > 1
                        lsubq = cqids.data{cqSub}(lqidx);
                        lsubq = lsubq{1};
                        if length(qidx) < lsubq
                            qidx = qidx(1);
                        else
                            qidx = qidx(lsubq);
                        end
                    end
                    qenum = stmData.data{ecol}(qidx);
                    if isnumeric(qenum)
                        qdata = enum2data(qenum);
                    end
                end % if length(qcidx)
                ocol = outCols.(make_valid_struct_key(qi));
                qtype = outData.datatype{ocol};
                qdata = sanitize_cell_value(qdata,qtype);
                if qtype == 's'
                    outData.data{ocol}{outRow} = qdata;
                else
                    outData.data{ocol}(outRow) = qdata;
                end
            end % for iq=1:length(qnums)
            cqids.data{cqBit} = cell(length(cqids.data{cqScq}),1);
          end % for in=1:nUt
        end % for istim = 1:nstim
    end % for isub = 1:nsub
else
    outData.data = sub_st.data;
    outRow = length(outData.data{1});
end

% % % % CONVERT outData.vars to meaningful codes
% 
if isfield(params.export,'var_name_map')
    vnmap = params.export.var_name_map;
    if isstruct(vnmap)
        % save original var names in outData.origvars
        outData.origvars = outData.vars;
        % convert var names
        for iv = 1:length(outData.vars)
            if isfield(vnmap,outData.vars{iv})
                outData.vars{iv} = vnmap.(outData.vars{iv});
            end
        end
    else
        warning(['vnmap is not a valid struct, therefore variables '...
            'will not be re-named']);
    end
end

% % print data out
% init file
fid = ensemble_init_fid(params);
fprintf(fid,'%s\n',cell2str(outData.vars,','));
if fid ~= 1
    % non-numeric cell delimiter
    if isfield(params.export,'non_num_delim')
        nndlm = params.export.non_num_delim;
    else
        nndlm = '';
    end
    % delimiter
    if isfield(params.export,'delimiter')
        dlm = params.export.delimiter;
    else
        % \t is the default delimiter
        dlm = '\t';
        params.export.delimiter = '\t';
    end
    for irow = 1:outRow
        row = '';
        for ivar = 1:length(outData.vars)
            item = outData.data{ivar}(irow);
            if iscell(item)
                if isempty(item)
                    item = '.';
                else
                    item = item{1};
                    if iscell(item)
                        item = cell2str(item);
                        % does this break anything??? FIXME
                        if isempty(item)
                            item = '.';
                        end
                    end
                end
            end
            if isnan(item)
                item = '.';
            end
            if isnumeric(item)
                item = num2str(item);
            end
            if islogical(item)
                if item
                    item = '1';
                else
                    item = '0';
                end
            end
            try item = regexprep(item,'\s','');
            catch
                error('error trying to remove white space from a cell')
            end
            item = regexprep(item,'[^a-zA-Z0-9.]','_');
            if ~isempty(regexp(item,'[^0-9.]','once')) && ~isempty(nndlm)
                item = [nndlm item nndlm];
            end
            if isempty(row)
                row = item;
            else
                row = sprintf('%s%s%s',row,dlm,item);
            end
        end
        fprintf(fid,'%s\n',row);
    end

    fclose(fid);

    % % % % GENERATE SAS IMPORT CODE
    % 
    if isfield(params,'sas') && isfield(params.sas,'fname')
        % expects params.sas.fname & params.fname
        fname = params.fname;
        sfname = params.sas.fname;
        sid = fopen(sfname,'wt');
        if sid == -1, error('Could not open SAS file: %s\n', sfname), end

        % Generate some header information
        fprintf(sid,'/* File generated by ensemble_export_respnstim.m, %s  */\n',datestr(now));
        if isfield(params.sas,'libpath')
            fprintf(sid,'LIBNAME xportlib ''%s'';\n',params.sas.libpath);
        end
        if isfield(params.sas,'options')
            fprintf(sid,'OPTIONS %s;\n',params.sas.options);
        end
        fprintf(sid,'\nDATA exported;\n')
        fprintf(sid,'  INFILE ''%s'' DLM=''%s'' FIRSTOBS=2;\n',fname,params.export.delimiter);
        %% FIXME: need some kind of logic to identify non-numeric columns and
        %% include a $ in the SAS INPUT command for that variable
        fprintf(sid,'  INPUT %s;\n\n',generate_sas_input_str(outData));
        fprintf(sid,'PROC CONTENTS;\nPROC MEANS mean stddev min max n;\nrun;\n');
        if isfield(params.sas,'libpath') && isfield(params.sas,'libsave') ...
                && isfield(params.sas,'libname')
            fprintf(sid,'DATA xportlib.%s;\nSET exported;\n',params.sas.libname);
        end
        fclose(sid);
    end
end

% % % % functions
% % sanitize_cell_value

function cval = sanitize_cell_value(cval,ctype)

if isempty(cval)
    if ctype ~= 's'
        cval = NaN;
    else
        cval = '.';
    end
end
if iscell(cval) && length(cval) == 1
    cval = cval{1};
end
if ctype ~= 's' && cval == '.'
    cval = NaN;
end

% % generate_sas_input_str

function input_str = generate_sas_input_str(struct)

input_str = '';
tempstr = '';

for iv=1:length(struct.vars)
    tempstr = sprintf('%s %s',tempstr,struct.vars{iv});
    if struct.datatype{iv} ~= 'n' && struct.datatype{iv} ~= 'l'
        tempstr = sprintf('%s $',tempstr);
    end
    if length(tempstr) > 100
        input_str = sprintf('%s%s\n',input_str,tempstr);
        tempstr='';
    end
end

if ~isempty(tempstr)
    input_str = sprintf('%s%s\n',input_str,tempstr);
end

input_str = sprintf('%s',input_str);
