function fsf=create_fsf
% fsf=create_fsf;
%
% Populates FEAT fsf structure with default values
%
% See also write_fsf

% 04/22/05 PJ


% FEAT version number
fsf.version=5.4;

% Analysis level
% 1 : First-level analysis
% 2 : Higher-level analysis
fsf.level=1;

% Which stages to run
% 0 : No first-level analysis (registration and/or group stats only)
% 7 : Full first-level analysis
% 1 : Pre-Stats
% 3 : Pre-Stats + Stats
% 2 :             Stats
% 6 :             Stats + Contrasts, Thresholding, Rendering
% 4 :                     Contrasts, Thresholding, Rendering
fsf.analysis=7;

% Delay before starting (hours)
fsf.delay=0;

% Use relative filenames
fsf.relative_yn=0;

% Balloon help
fsf.help_yn=1;

% Run Featwatcher
fsf.featwatcher_yn=1;

% Cleanup first-level standard-space images
fsf.sscleanup_yn=1;

% Output directory
fsf.outputdir='';

% TR(s)
fsf.tr=[];

% Total volumes
fsf.npts=[];

% Delete volumes
fsf.ndelete=0;

% Number of first-level analyses
fsf.multiple=1;  % Not a flag - actual num of analyses

% Higher-level input type
% 1 : Inputs are lower-level FEAT directories
% 2 : Inputs are cope images from FEAT directories
fsf.inputtype=1;

% Carry out pre-stats processing?
fsf.filtering_yn=1;

% Brain/background threshold, %
fsf.brain_thresh=10;

% Post-stats-only directory copying
% 0 : Overwrite original post-stats results
% 1 : Copy original FEAT directory for new Contrasts, Thresholding, Rendering
fsf.newdir_yn=0;

% Slice timing correction
% 0 : None
% 1 : Regular up (0, 1, 2, 3, ...)
% 2 : Regular down
% 3 : Use slice order file
% 4 : Use slice timings file
% 5 : Interleaved (0, 2, 4 ... 1, 3, 5 ... )
fsf.st=0;

% Slice timings file
fsf.st_file='';

% Motion correction
% 0 : None
% 1 : MCFLIRT
fsf.mc=1;

% Spin-history (currently obsolete)
fsf.sh_yn=0;

% BET brain extraction
fsf.bet_yn=1;

% Spatial smoothing FWHM (mm)
fsf.smooth=5;

% Intensity normalization
fsf.norm_yn=0;

% Highpass temporal filtering
fsf.temphp_yn=0;

% Lowpass temporal filtering
fsf.templp_yn=0;

% MELODIC ICA data exploration
fsf.melodic_yn=0;

% Carry out main stats?
fsf.stats_yn=1;

% Carry out prewhitening?
fsf.prewhiten_yn=1;

% Higher-level modelling
% 3 : Fixed effects
% 0 : Mixed Effects: Simple OLS
% 2 : Mixed Effects: FLAME (stage 1 only)
% 1 : Mixed Effects: FLAME (full)
fsf.mixed_yn=1;

% Number of EVs
fsf.evs_orig=0;
fsf.evs_real=0;

% Number of contrasts
fsf.ncon_orig=0;
fsf.ncon_real=0;

% Number of F-tests
fsf.nftests_orig=0;
fsf.nftests_real=0;

% Add constant column to design matrix? (obsolete)
fsf.constcol=0;

% Carry out post-stats steps?
fsf.poststats_yn=1;

% Pre-threshold masking?
fsf.threshmask='';

% Thresholding
% 0 : None
% 1 : Uncorrected
% 2 : Voxel
% 3 : Cluster
fsf.thresh=3;

% P threshold
fsf.prob_thresh=0.05;

% Z threshold
fsf.z_thresh=2.3;

% Z min/max for colour rendering
% 0 : Use actual Z min/max
% 1 : Use preset Z min/max
fsf.zdisplay=0;

% Z min in colour rendering
fsf.zmin=1;

% Z max in colour rendering
fsf.zmax=15;

% Colour rendering type
% 0 : Solid blobs
% 1 : Transparent blobs
fsf.rendertype=1;

% Background image for higher-level stats overlays
% 1 : Mean highres
% 2 : First highres
% 3 : Mean functional
% 4 : First functional
% 5 : Standard space template
fsf.bgimage=1;

% Registration?
fsf.reg_yn=1;

% B0 fieldmap unwarping?
fsf.regunwarp_yn=0;

% Dwell/Asymmetry ratio 
fsf.dwellasym=-1;

% Registration to initial structural
fsf.reginitial_highres_yn=0;

% Search space for registration to initial structural
% 0   : No search
% 90  : Normal search
% 180 : Full search
fsf.reginitial_highres_search=90;

% Degrees of Freedom for registration to initial structural
fsf.reginitial_highres_dof=12;

% Do nonlinear registration to initial structural?
fsf.reginitial_highres_nonlinear_yn=0;

% Registration to main structural
fsf.reghighres_yn=0;

% Search space for registration to main structural
% 0   : No search
% 90  : Normal search
% 180 : Full search
fsf.reghighres_search=90;

% Degrees of Freedom for registration to main structural
fsf.reghighres_dof=12;

% Do nonlinear registration to main structural?
fsf.reghighres_nonlinear_yn=0;

% Registration to standard image?
fsf.regstandard_yn=1;

% Standard image
fsf.regstandard='/usr/local/fsl/etc/standard/avg152T1_brain.hdr';

% Search space for registration to standard space
% 0   : No search
% 90  : Normal search
% 180 : Full search
fsf.regstandard_search=90;

% Degrees of Freedom for registration to standard space
fsf.regstandard_dof=12;

% Do nonlinear registration to standard space?
fsf.regstandard_nonlinear_yn=0;

% High pass filter cutoff
fsf.paradigm_hp=-1;

% Number of lower-level copes feeding into higher-level analysis
fsf.ncopeinputs=0;

% 4D AVW data or FEAT directory (1)
fsf.feat_files={};

% Insert EV structure creation here
fsf.ev = create_ev;

fsf.con = [];  % create contrasts with create_con
fsf.Ftest = [];  %

% Contrast & F-tests mode
% real : control real EVs
% orig : control original EVs
fsf.con_mode_old='orig';
fsf.con_mode='orig';

% Read in default contrast structures if desired

% Contrast masking - use >0 instead of thresholding?
fsf.conmask_zerothresh_yn=0;

% This is an EV X EV matrix of what contrast should be masked with what.
% Mask contrast index in row with contrast index in column
% e.g. Mask real contrast/F-test 1 with real contrast/F-test 2?
fsf.conmask_mtx = [];

% Do contrast masking at all?
fsf.conmask1_1=0;
