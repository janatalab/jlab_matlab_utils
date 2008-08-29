function [ave_chdr] = sesfhdr2avetspecs(ses_fhdr, ave_chdr)
%function [ave_chdr] = sesfhdr2avetspecs(ses_fhdr, ave_chdr)
%
%Copies run info, etc to tspecs in ave file
%

% Modification history
% 7/20/95 -- PJ 	Created routine
% 7/22/95 -- PJ		added tests for size of TSpec in average chdr
%			and value in LSpec
%			The LSpec must be at least 28 in order for wt_ave_hdr_v to succeed

ses_hdr_offsets_v; ave_hdr_offsets_v;

if size(ave_chdr,2) < 19 
  ave_chdr = ave_chdr + zeros(size(ave_chdr,1), 19 - size(ave_chdr,2));
end

if any(ave_chdr(:,LSpec) < 28)
  ave_chdr(:,LSpec) = ones(size(ave_chdr,1),1) * 28;
end
nrows = size(ave_chdr,1);
ave_chdr(:,RunDateMo_Ave) = ones(nrows,1)*ses_fhdr(RunDateMo);
ave_chdr(:,RunDateDay_Ave) = ones(nrows,1)*ses_fhdr(RunDateDay);
ave_chdr(:,RunDateYr_Ave) = ones(nrows,1)*ses_fhdr(RunDateYr);
ave_chdr(:,RunTimeHr_Ave) = ones(nrows,1)*ses_fhdr(RunTimeHr);
ave_chdr(:,RunTimeMin_Ave) = ones(nrows,1)*ses_fhdr(RunTimeMin);
ave_chdr(:,RunTimeSec_Ave) = ones(nrows,1)*ses_fhdr(RunTimeSec);
ave_chdr(:,SubjID_Ave) = ones(nrows,1)*ses_fhdr(SubjID);
ave_chdr(:,Handed_Ave) = ones(nrows,1)*ses_fhdr(Handed);
ave_chdr(:,Sex_Ave) = ones(nrows,1)*ses_fhdr(Sex);
ave_chdr(:,Age_Ave) = ones(nrows,1)*ses_fhdr(Age);
ave_chdr(:,ExperID_Ave) = ones(nrows,1)*ses_fhdr(ExperID);
ave_chdr(:,EdVer_Ave) = ones(nrows,1)*ses_fhdr(EdVer);
ave_chdr(:,CalFlag_Ave) = ones(nrows,1)*ses_fhdr(CalFlag);
