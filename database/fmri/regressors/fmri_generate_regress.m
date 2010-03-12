function [out] = fmri_generate_regress(pinfo, minfo, sess)
% generates regressors for SPM5 or FSL fMRI analysis design matrices
% 
%   [out] = fmri_generate_regress(pinfo,minfo,sess)
% 
% Given a set of parameters, this function generates regressors by calling
% sub-functions named "fmri_regress_%regid%.m", where 'regid' is the name
% of the regressor as provided in the model information structure. Upon
% compiling a set of regressors, this function will package those
% regressors as necessary, for either FSL or SPM5, and will return the
% requisite structure (an fsf.ev structure for FSL, or a .cond or .regress
% part of the .sess structure for SPM).
% 
% REQUIRES
%   pinfo - run information structure, including presentation data
%   minfo - model information structure
%   sess  - session information structure
% 
% RETURNS
%   out - structure in either fsf.ev, spmsess.regress, or spmsess.cond
%       format
% 
% 2009.11.05 FB - adapted from fmri_fsl_generate_evs,
% fmri_spm_generate_regress, and fmri_spm_generate_conds

% generate final structure for SPM or FSL?
try USE_SPM = pinfo.USE_SPM; catch USE_SPM = 0; end
try USE_FSL = pinfo.USE_FSL; catch USE_FSL = 0; end
if (~USE_SPM && ~USE_FSL) || (USE_SPM && USE_FSL)
  error('must choose either FSL or SPM for output');
end

% compile condition
condlist = fieldnames(minfo.condlist);
ncond = length(condlist); % number of conditions

% iterate over conditions, form condition list
cnames = {};
consets = {};
durations = {};

for icond=1:ncond

  cond = condlist{icond};

  cstr = sprintf('fmri_cond_%s',cond);
  if ~exist(cstr,'file')
    warning('Unknown regressor: %s',cond);
    continue
  end

  cfun = str2func(cstr);

  [lons,ldurs] = cfun(pinfo,minfo,sess);

  if isempty(lons) || isempty(lons) || (length(lons) ~= length(ldurs))
    warning('no regressors generated for regressor %s',cond);
    continue
  end

  cnames = [cnames cond];
  consets = [consets lons];
  durations = [durations ldurs];

end

% compile parametric condition modulations
pmodlist = minfo.pmodlist;
nmod = size(pmodlist,1); % number of parametric modulation sets
pmods = cell(length(cnames),1);
ccm = parse_fh(minfo.cond_cue_map);
pc = set_var_col_const(pinfo.vars);

% iterate over parametric modulations
for imod = 1:nmod
  cond_name = pmodlist{imod,1};
  
  % Determine which conds array element we are attaching the parametric
  % modulation to.
  cond_idx = strmatch(cond_name,cnames,'exact');
  if ~isempty(cond_idx)
    modlist = pmodlist{imod,2};
    nlmod = length(modlist);
    for ilmod = 1:nlmod
      modname = modlist{ilmod};
      lp = struct();
      lp.name = modname;

      cue_type = ccm(modname);

      rfilt.include.all = cue_type.filt;
      rdata = ensemble_filter(pinfo,rfilt);
	  
      % Get the parameter value for each event
      if isempty(rdata.data{1})
  	    fprintf('Could not match desired parametric modulator: %s\n', modname);
      else
        % Extract the button codes
        respvect = rdata.data{pc.RESP_CODE};

        % Map the button codes to the desired values
        nval = length(respvect);
        % Get correct index in response mappings
        rmap = pinfo.resp_mapping;
        mapidx = strmatch(cue_type.name, {minfo.resp_mappings{rmap}{:,1}},'exact');
        param = [];
        for ival = 1:nval
          currmask = minfo.resp_mappings{rmap}{mapidx,2} == str2num(respvect{ival}{1});
          param(ival) = minfo.resp_mappings{rmap}{mapidx,3}(currmask);
        end
	    lp.param = param(:);
        pmods{cond_idx} = [pmods{cond_idx} lp];
      end % if ~isempty(ridx);
    end % for imod
  end % if ~isempty(cond_idx)
end % for icondimod

% compile regressors
reglist = minfo.reglist;
nreg = length(reglist);  % Number of regressors
reg_template = struct('name','','val',[]);

names = {};
vals = [];

% iterate over regressors, form regressor list
for ireg = 1:nreg
    
  regid = reglist{ireg};

  linfo = pinfo;
  linfo.regid = regid;
  
  if strmatch('stim_',regid)
    [lnames,lvals] = fmri_regress_stim(linfo,minfo,sess);
  else
    regstr = sprintf('fmri_regress_%s',regid);
    if ~exist(regstr,'file')
      warning('Unknown regressor: %s',regid);
      continue
    end
    
    regfun = str2func(regstr);
  
    [lnames,lvals] = regfun(linfo,minfo,sess);
  end % if strmatch('stim_
  
  if isempty(lnames) || isempty(lvals)
    warning('no regressors generated for regressor %s',regid);
    continue
  end
  
  names = [names lnames];
  
  vrows = size(vals,1);
  lrows = size(lvals,1);
  if vrows < lrows && vrows > 0
    vals(vrows+1:lrows,:) = 0;
  elseif lrows < vrows
    lvals(lrows+1:vrows,:) = 0;
  end % if vrows < lrows
  
  vals = [vals lvals];
  
end % for ireg=1:nreg

% % % BUILD TOOL-SPECIFIC REGRESSOR STRUCTURE

if USE_FSL
    
  evoutdir = pinfo.evoutdir;
  evstub = pinfo.evstub;

  ev = create_ev;
  nev = 0;
  % add regressors
  for ir=1:length(names)
    nev = nev+1;
    outfname = fullfile(evoutdir,sprintf(evstub,nev));
    ev(nev).name = names{ir};
    ev(nev).fname = outfname;
    ev(nev).shape = 2;
    regdata = vals{ir};
    save(outfname,'regdata','-ascii');
  end
  
  % add conditions
  ortho_idxs = [];
  for ic=1:length(cnames)
    nev = nev+1;
    ev(nev).name = cnames{ic};
    ev(nev).shape = 0;
    ev(nev).convolve = 3;
    ev(nev).on = consets{ic};
    ev(nev).off = consets{ic} + durations{ic};
    if isempty(pmods{ic}), continue, end
    lpmods = pmods{ic};
    lortho_idxs = [];
    for ip=1:length(lpmods)
      % add parametric modulations
      nev = nev+1;
      outfname = fullfile(evoutdir,sprintf(evstub,nev));
      ev(nev).name = lpmods{ip}.name;
      ev(nev).fname = outfname;
      ev(nev).shape = 2;
      regdata = build_regress(consets{ic},durations{ic},...
          lpmods{ip}.param,pinfo.TR,pinfo.dt,pinfo.actual_nvol);
      save(outfname,'regdata','-ascii');
      lortho_idxs = [ortho_idxs nev];
    end
    ortho_idxs = [ortho_idxs lortho_idxs];
  end

  [ev.ortho] = deal(zeros(1,nev));
  
  % set orthogonalization for pmod regressor sets
  for io=1:length(ortho_idxs)
    lortho = ortho_idxs{io};
    for il=1:length(lortho)
      ev(lortho(il)).ortho(lortho) = 1;
    end
  end
  
  out = ev;
  
elseif USE_SPM
    
  % Initiate the condition output variable
  if ~isempty(consets)
    condstruct = struct('name','','onset',[],'duration',[],'pmod',[]);
    for ic=1:length(cnames)
      condstruct(ic).name = cnames{ic};
      condstruct(ic).onset = consets{ic};
      condstruct(ic).duration = durations{ic};
      condstruct(ic).pmod = pmods{ic};
    end
    [condstruct(1:end).tmod] = deal(0);
    out.cond = condstruct;
  end
  
  % Initialize the regressors output variable to an empty instance of the structure
  if ~isempty(vals)
    regressors = struct('name','','val',[]);
    for ir=1:length(names)
      regressors(ir).name = names{ir};
      regressors(ir).val = vals(:,ir);
    end
    out.regress = regressors;
  end
  
end % if USE_FSL/USE_SPM
