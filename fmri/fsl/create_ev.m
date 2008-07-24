function ev=create_ev
% ev=create_ev;
%
% Creates a regressor for use in FEAT model

% 2005 Petr Janata


% Basic waveform shape (EV 1)
% 0 : Square
% 1 : Sinusoid
% 2 : Custom (1 entry per volume)
% 3 : Custom (3 column format)
% 4 : Interaction
ev.shape=0;

% Convolution (EV 1)
% 0 : None
% 1 : Gaussian
% 2 : Gamma
% 3 : Double-Gamma HRF
% 4 : Gamma basis functions
% 5 : Sine basis functions
% 6 : FIR basis functions
ev.convolve=0;

% Convolve phase (EV 1)
ev.convolve_phase=0;

% Apply temporal filtering (EV 1)
ev.tempfilt_yn=0;

% Add temporal derivative (EV 1)
ev.deriv_yn=0;

%
% Different EV shapes take different sets of parameters
%

%
% Parameters for square waves
% 

% Skip (EV 1)
ev.skip=0;

% Off (EV 1)
ev.off=[];

% On (EV 1)
ev.on=[];

% Phase (EV 1)
ev.phase=[];

% Stop (EV 1)
ev.stop=[];

%
% Parameters for custom EVs
%
ev.fname='';

% The orthogonalization vector has to obviously be set once the EVs are all known
% Orthogonalise EV 1 wrt EV 0
ev.ortho = [];

%
% Convolution parameters
%
ev.gausssigma = 2.8;
ev.gaussdelay = 3;