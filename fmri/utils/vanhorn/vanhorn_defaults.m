% vanhorn_defaults.m
%
% Default variable definitions used in Jack Van Horn's analysis scripts
%


% ***************************************
% Define Endian Type for Image Data Input
ENDIANIN='ieee-le';   % Use this when running on Intel-Based Processors and
                      % with data written using an Intel Processor
                      % e.g. A DELL Linux Box

%ENDIANIN='ieee-be';  % Use this when running on SUN-Based Processors and
                      % with data written using a SUN Processor
                      % e.g. dexter.dartmouth.edu

% ****************************************
% Define Endian Type for Image Data Output
ENDIANOUT='ieee-le';  % Use this when running on Intel-Based Processors and
                      % with data written using an Intel Processor
                      % e.g. A DELL Linux Box

%ENDIANOUT='ieee-be'; % Use this when running on SUN-Based Processors and
                      % with data written using a SUN Processor
                      % e.g. dexter.dartmouth.edu

		      
% *********************
% Set Masking Threshold
THRESH=550.0;

% **************************************
% Flag for Normalizing the Design Matrix
NORMDESMTX=0;

% ***************************
% fMR Data Normalization Mode
DATA_NORM_MODE = 3;

% ***************************
% Flag for Outlier Protection
OUTLIER_CORRECTION=1;
OUTLIER_THRESH=4.0;

% ************************************************
% Flag for The Existence of a Matlab Contrast File
CONTRAST_FILE = 0;

% ************************************
% Critical Probability Threshold Value
pcrit=0.05;

% ***************************
% Set A Small Tolerance Value
TINY=1e-5;

% ********************
% F_indiv Scatterplots
F_INDIV_PLOTS=0;

% ***************
% Local Max Plots
LOCAL_MAX_PLOTS=1;

% ************************************************************
% If Christian Buchel's Structural Equation Modeling Routine's
% are Available set SEM=1, else SEM=0
SEM=0;
