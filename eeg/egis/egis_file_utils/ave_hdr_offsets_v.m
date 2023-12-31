%ave_hdr_offsets.m
%
%Declares a set of global constants for writing the fields of the
%EGIS average header vectors.
%
%Based on ses_hdr_offsets.m
%
%See also:rd_ses_hdr_vecs, wt_ses_hdr_vecs, wt_ave_hdr_v.

%Modification History:
%	06/06/95 PJ Created file

% The fields given below reference positions in the fhdr information vector

%Average File Description from "EGIS Hdr Conv 1.02"
BytOrd=1;
HdrVer=2;
LHeader=3;
LData=4;

%Average File Session Information
% Fields 5 through 10 are used for storing run date and run time information
% in session files.  They are not used in average files
LastDone=11;  	% Number of subjects in average file
ScaleBins=12;		
ScaleCal=13;
BaseDur=14;
% Fields 15 through 17 are empty in average files

%Average File Data Description
NCells=18;
NChan=19;
LComment=20;
LText=21;
LPad=22;
BrdGain=23;
LCellHdr=24;

% The fields given below reference positions in the chdr information matrix

%Average Cell Offsets
CellID=1;
NSubjects=2;
NObs=2;
NPoints=3;
NSamps=3;
NSamp=3;
SampRate=4;
LSpec=5;

%Subject specifics array
NAvg=6;		% number of trials in cell average
RunDate_Ave=7;
RunDateMo_Ave=7;
RunDateDay_Ave=8;
RunDateYr_Ave=9;
RunTime_Ave=10;
RunTimeHr_Ave=10;
RunTimeMin_Ave=11;
RunTimeSec_Ave=12;
SubjID_Ave=13;
Handed_Ave=14;
Sex_Ave=15;
Age_Ave=16;
ExperID_Ave=17;
EdVer_Ave=18;
CalFlag_Ave=19;

