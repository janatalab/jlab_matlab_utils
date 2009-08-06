function fsf = set_fsf_defaults(fsf, setname)
% fsf = set_fsf_defaults(fsf);
%
% Sets fields of the FEAT analysis structure (fsf) according to the defaults
% associate with a particular set of settings.
%
% setname - name of set 
%

% 11/03/05 Petr Janata

try fsf; catch fsf=[]; end
try setname; catch setname = ''; end

if isempty(setname) | isempty(fsf)
  fsf = create_fsf;
end

switch setname
  case {'remove_physio','remove_variance'}
    fsf.analysis = 2;  % 2=stats only
    fsf.poststats_yn = 0;
    fsf.paradigm_hp = 0;
    fsf.prewhiten_yn = 0; % no prewhitening
    fsf.filtering_yn = 0; % no filtering
    fsf.mc = 0;  % turn motion-correction off
    fsf.bet_yn = 0; % no brain extraction
    fsf.smooth = 0; % no smoothing
    fsf.regstandard_yn = 0;  % no registration
    
    fsf.con_mode_old = 'real';
    fsf.con_mode = 'real';
    
  case 'evaluate_physio'
    fsf.paradigm_hp = 0;
    fsf.prewhiten_yn = 0; % no prewhitening
    fsf.filtering_yn = 0; % no filtering
    fsf.analysis = 6;  % 4=Stats + Contrasts, Thresholding, Rendering
    fsf.thresh=2;  % 0=None, 1=Uncorrected, 2=Voxel, 3=Cluster
    fsf.regstandard_yn = 0;  % no registration
    fsf.mc = 0;
    fsf.bet_yn = 0;
    fsf.smooth = 0;
    fsf.con_mode_old = 'real';
    fsf.con_mode = 'real';
    fsf.stats_yn = 1;
    fsf.poststats_yn = 1;
    fsf.featwatcher_yn = 1;
    fsf.tsplot_yn = 1;
    
  case 'coreg'
    fsf.analysis = 0;  % Registration
    fsf.reginitial_highres_yn=1; % coplanar
    fsf.reghighres_yn=1; % high res
    fsf.regstandard_yn = 1;
    fsf.stats_yn = 0;
    fsf.poststats_yn = 0;
    fsf.reg_yn=1;
    
  case 'subject_sess_run'
    fsf.level = 2;  % 2=higher level analysis
    fsf.analysis = 6; % 6=stats+post-stats
    fsf.inputtype = 1;  % 1 : Inputs are lower-level FEAT directories
end
