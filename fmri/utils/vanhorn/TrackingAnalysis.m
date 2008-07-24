% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Script for Performing Slice-Based
% Multiple Regression on fMRI Time Series
%
% Written by:
% John Darrell Van Horn, Ph.D., M.Eng.
% Dartmouth Brain Imaging Center
% Dartmouth College
% Hanover, NH 03753
%
% Phone: (603) 646-2909
% Email: John.D.Van.Horn@dartmouth.edu
%
%
% Last modified: 11-29-2001
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% **********************
% Clear Matlab Workspace
clear;
clc;

% ****************
% Close All Figure
close all;

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

% *************************
% Start Matlab Elapsed Time
tic;

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

% *********************
% Image File Dimensions
TR=2.0;
RUNS=4;
VOLS_PER_RUN=126;
XDIM=64;
YDIM=64;
ZDIM=25;
XVOXDIM=3.75;
YVOXDIM=3.75;
ZVOXDIM=5.5;
N=RUNS*VOLS_PER_RUN;
NUMVOX=XDIM*YDIM*ZDIM;

% ************************************************************
% If Christian Buchel's Structural Equation Modeling Routine's
% are Available set SEM=1, else SEM=0
SEM=1;

% *************************
% Select Design Matrix File
[DesMtxFileName, CWD]=uigetfile('*.mat','Select Design Matrix Matlab File');

% *********************
% Load Design Matrix, X
disp('Loading Design Matrix...');
close;
eval(['load ' [CWD DesMtxFileName] ';']);
DesMtx=DesMtxB;  
clear DesMtxA DesMtxB;
D=size(DesMtx);  %For the rest of the script the Design Matrix is called DesMtx
VARIABLES=D(1);
TIMEPOINTS=D(2);
if(length(D)==3)
  SLICES=D(3);
elseif(length(D)==2)
  SLICES=1;
else
  disp(['This ain`t no stinkin design matrix!']);
  exit;
end

% ********************************************************************
% If Design Matrix is not 3D Slice-based (i.e. is only a single slice)
% Replicate it to create a Design volume
if(SLICES==1)
     tmp=DesMtx;
     DesMtx=zeros(VARIABLES,TIMEPOINTS,ZDIM);
     for i=1:ZDIM,
	 DesMtx(:,:,i)=tmp(:,:);
     end
     D=size(DesMtx);
     VARIABLES=D(1);
     TIMEPOINTS=D(2);
     SLICES=ZDIM;
end

% *************************************************
% Create Selection Vector and Define Masking Vector 
% for Removing Effects of Confounding Variables
% e.g. those variables that will be regressed out 
% for the reduced model
S=zeros(1,VARIABLES);
S([37:end])=1;

% **********************
% Load a Contrast Matrix
if(CONTRAST_FILE==1),

  % Select Contrast Matrix Matlab File
  [ContMtxFileName, CWD]=uigetfile('*.mat','Select Contrast Matrix Matlab File');
  eval(['if(exist(' [CWD ContMtxFileName] '))']);
  eval(['load ' [CWD ContMtxFileName] ';']);
else
  disp(['No contrast matrix exists by that name.']);
  disp(['Carrying on anyway......']);
  C=[(1-eye(VARIABLES-sum(S),VARIABLES-sum(S))) ones(VARIABLES-sum(S),sum(S))];
end

% ***************************
% Select Image List Text File
[ImageFileList, DataDir]=uigetfile('*.txt');

% ***********************
% Read in Image Filenames
disp('Reading in Filenames...');
close;
drawnow;
fname=[DataDir ImageFileList];
fid=fopen(fname,'r',ENDIANIN);

for i=1:N,
  F=fgetl(fid);
  FILES(i,[1:length(F)])=F;
end
fclose(fid);

% *************************************
% Read in Functional Image Data Volumes
disp('Loading Functional Image Data...');
Y=zeros(N,NUMVOX);
for i=1:N,
	fid=fopen(FILES(i,:),'r',ENDIANIN);
        U=fread(fid,NUMVOX,'int16');
        Y(i,:)=U';
        fclose(fid);
end

% *************************
% Load or Create Image Mask
if(exist('Mask.img'))
  disp('Loading mask image from disk...');
  fid=fopen('Mask.img','r',ENDIANIN);
  mask=fread(fid,NUMVOX,'int8');
  fclose(fid);
else
  disp('Generating Image Mask...');
  I=(Y>THRESH);
  mask=prod(I);
  clear I;

  % Save Mask Image to PWD
  fid=fopen(['Mask.img'],'w',ENDIANOUT);
  fwrite(fid,mask,'int8');
  fclose(fid);
  str=['!makeheader4 Mask.hdr ', num2str(XDIM), ' ', num2str(YDIM), ' ', num2str(ZDIM),...
  ' ', num2str(XVOXDIM), ' ', num2str(YVOXDIM), ' ', num2str(ZVOXDIM), ' 1 8 1 0 1 0 2'];
  eval(str);
end
V=find(mask>0);

if(length(V)>0)
  disp(['Analysing ' num2str(sum(length(V))) ' voxels above threshold...']);
else
  disp(['No voxels found above threshold.']);
  clear;
  quit;
end

% ********************************
% Standard Normalize Design Matrix
if (NORMDESMTX==1)
   disp('Normalizing Design Matrix...');
   Z=stdnormalize(DesMtx);
else
   Z=DesMtx;
end

% ******************************
% Normalize fMR Image Data by... 
switch DATA_NORM_MODE
  case 1     % 1) Proportional Scaling

     disp('Normalizing Image Data by Proportional Scaling...');
     ImageMeans=zeros(1,N);
     ImageMeans=mean(Y(:,V)');
     for i=1:N,
	for j=1:length(V),
            Y(i,V(j))=Y(i,V(j))./ImageMeans(i);
        end
     end

  case 2     % 2) Standard Normalization Over All Scans 
     disp('Normalizing Image Data by Standard Normalization...');
     ImageMeans=zeros(1,N);
     ImageMeans=mean(Y(:,V)');
     ImageStds=std(Y(:,V)');
     for i=1:N,
        for j=1:length(V),
            Y(i,V(j))=(Y(i,V(j))-ImageMeans(i))./ImageStds(i);
        end
     end

  case 3     % 3) Standard Normalization of Scans Within Run 
     disp('Normalizing Image Data by Standard Normalization Within Run...');

     for i=1:RUNS,
     ImageMeans=zeros(1,N);
     ImageMeans=mean(Y([(i-1)*VOLS_PER_RUN+1:i*VOLS_PER_RUN],V)');
     ImageStds=std(Y([(i-1)*VOLS_PER_RUN+1:i*VOLS_PER_RUN],V)');
        for k=((i-1)*VOLS_PER_RUN+1):(i*VOLS_PER_RUN),
           for l=1:length(V),
            Y(k,V(l))=(Y(k,V(l))-ImageMeans(i))./ImageStds(i);
           end
        end
     end

  otherwise
     disp(['Unknown MR Data Normalization Method']);
     exit;
end

     save Y;

% **************************************
% Perform Slice-Wise Pseudo-Inversion of 
% the Design Matrix for Use Below
PinvDesMtx = zeros(VARIABLES,VARIABLES,SLICES);
for i=1:SLICES,
    PinvDesMtx(:,:,i)=pinv(Z(:,:,i)*Z(:,:,i)');
end


% **************************************
% Perform Slice-wise Multiple Regression
disp('Performing Slice-wise Multiple Regression...');
Beta=zeros(VARIABLES,NUMVOX);
H=zeros(1,NUMVOX);
E=zeros(1,NUMVOX);
T=zeros(1,NUMVOX);
Lambda=zeros(1,NUMVOX);
Lambda_reduced=zeros(1,NUMVOX);
F_omnibus=zeros(1,NUMVOX);
F_reduced=zeros(1,NUMVOX);
F_reduced_mask=zeros(1,NUMVOX);
F_indiv=zeros(VARIABLES-sum(S),NUMVOX);

if(OUTLIER_CORRECTION==1)
     disp(' ');
     disp('Outlier Correction in effect...');
     disp(' ');
end

zlast=0;
for i=1:length(V),

  % When enabled, inspect y vector and set outliers equal to zero
  % and set values in the data equal to mean of the rest of the series
  if(OUTLIER_CORRECTION==1)  
    q=Y(:,V(i));
    q=(q-mean(q))./std(q);
    for j=1:N
      if(abs(q(j))>OUTLIER_THRESH)
         Y(j,V(i))=(sum(Y(:,V(i)))-Y(j,V(i)))./(N-1);
      end
    end
  end

  y=Y(:,V(i));
  y=y-mean(y);

  z=floor(V(i)./(XDIM*YDIM))+1;
  if (z~=zlast)
     x=Z(:,:,z);
     zlast=z;
     disp(['Now analyzing data on slice ' num2str(z) '...']);
  end

  b=PinvDesMtx(:,:,z)*x*y;
  Beta(:,V(i))=b;
  yprime_full=x'*Beta(:,V(i));
  T(V(i))=(y-mean(y))'*(y-mean(y));
  E(V(i))=(y-yprime_full)'*(y-yprime_full);
  H(V(i))=T(V(i))-E(V(i));

  Lambda(V(i))=det(E(V(i)))./det(T(V(i)));
  if((Lambda(V(i))<0.0)|(Lambda(V(i))>1.0))
     disp(['Lambda out of range on slice ' num2str(z) '.']);
     Lambda(V(i))=1.0;
  end

  F_omnibus(V(i))=((N-VARIABLES-1)./VARIABLES).*((1-Lambda(V(i)))./Lambda(V(i)));
  
  yprime_reduced=x'*(Beta(:,V(i)).*S');
  EplusH_reduced=(y-yprime_reduced)'*(y-yprime_reduced);
  Lambda_reduced(V(i))=det(E(V(i)))./det(EplusH_reduced);
  if((Lambda_reduced(V(i))<0.0)|(Lambda_reduced(V(i))>1.0))
     disp(['Lambda_reduced out of range on slice ' num2str(z) '.']);
     Lambda_reduced(V(i))=1.0;
  end
  F_reduced(V(i))=((N-VARIABLES-1)./(VARIABLES-sum(S))).*...
((1-Lambda_reduced(V(i)))./Lambda_reduced(V(i)));


[CONTRASTS A]=size(C);


  for j=1:CONTRASTS,
     if(sum(C(j,:))>0)                           % e.g. Sum of coeffs GT 0.0, assume a Reduced Model Vector
       yprime_j=x'*(Beta(:,V(i)).*C(j,:)');
       EplusH_j=(y-yprime_j)'*(y-yprime_j);
       Lambda_j=det(E(V(i)))./det(EplusH_j);
     elseif(abs(sum(C(j,:)))<TINY)               % e.g. Abs sum of coeffs <= TINY, Assume a Contrast Model Vector
       yprime_j=x'*(Beta(:,V(i)).*C(j,:)');
       H_j=yprime_j'*yprime_j;
       E_j=y'*y-yprime_j'*yprime_j;
       EplusH_j=(E_j+H_j);
       Lambda_j=det(E_j)./det(EplusH_j);
     else
       disp(['Contrast ' num2str(j) ' appears to be invalid.  Please check.']);
     end

     if((Lambda_j<0.0)|(Lambda_j>1.0))
        disp(['Lambda_j ' num2str(Lambda_j) ' out of range on slice ' num2str(z) '.']);
        Lambda_j=1.0;
     end
     
     F_indiv(j,V(i))=((N-VARIABLES-1)./1).*((1-Lambda_j)./Lambda_j);

     if(F_INDIV_PLOTS==1)
       if(Lambda_j<0.5)
         subplot(1,2,1);
         plot(x(j,:),(y-yprime_j),'.');
         hold on;
         y_reg=x(j,:).*Beta(j,V(i));
         plot(x(j,:),y_reg,'r');
         axis square;
         xlabel('X_j');
         ylabel('Y_i');
         title(['F_indiv #' num2str(j) ' = ' num2str(F_indiv(j,V(i))) ' (Lambda_j = ' num2str(Lambda_j) ') on slice ' num2str(z) ' at voxel ' num2str(V(i))]);
         hold off;
         subplot(1,2,2);
         plot((y-yprime_j),'g.');
         hold on;
         plot(y_reg,'r');
	 xlabel('Time');
         ylabel('Y_i');
         title(['Beta_j = ' num2str(Beta(j,V(i))) ', R-sqr = ' num2str(1-Lambda_j)]);
         hold off;
         drawnow;
         pause(0.5);
         print VoxelPlots.eps -dpsc -append;
      end
     end

  end

  % Periodically display program elapsed and projected time info to stdout
  if(mod(i,XDIM*YDIM)==0)
    pctcomplete = floor(100.*i./length(V));
    disp([ num2str(pctcomplete) ' percent complete...']);
    disp(['Program Elapsed Running Time = ' num2str(floor(toc/60)) ':' num2str(60*mod((toc/60),1)) ]);
    runtime=toc./(pctcomplete/100);
    timeleft=runtime-toc;
    disp(['Projected to run for ' num2str(floor(timeleft/60)) ':' num2str(60*mod((timeleft/60),1)) ' more minutes']);
  end

end

% **********************************
% Critical F Value for Reduced Model
Fcrit=Fcritical(pcrit,(VARIABLES-sum(S)),(N-VARIABLES-1));
Rsqr_crit=(Fcrit*(VARIABLES-sum(S))./(N-VARIABLES-1))./(1+(Fcrit*(VARIABLES-sum(S))./(N-VARIABLES-1)));
disp([' ']);
disp(['Critical Reduced Model F(' num2str(VARIABLES-sum(S)) ',' num2str(N-VARIABLES-1) ')=' num2str(Fcrit) ' at p=' num2str(pcrit)]);
disp(['Model R-Squared at Critical F = ' num2str(Rsqr_crit)]);
disp([' ']);

% Threshold F_Reduced Image at Critical F
%F_reduced_mask(V)=((F_reduced(V)==Fcrit)|(F_reduced(V)>Fcrit));
%for i=1:CONTRASTS,
%      F_indiv(i,:)=F_indiv(i,:).*F_reduced_mask;
%end

% *******************
% Save Results to PWD
disp('Saving Results Images...');

% ************************
% Beta Parameter Estimates
for i=1:VARIABLES,
fid=fopen(['Beta_' num2str(i) '.img'],'w',ENDIANOUT);
fwrite(fid,Beta(i,:),'float32');
fclose(fid);
str=['!makeheader4 ', ['Beta_' num2str(i) '.hdr '] , num2str(XDIM),...
' ', num2str(YDIM), ' ', num2str(ZDIM),...
' ', num2str(XVOXDIM), ' ', num2str(YVOXDIM), ' ', num2str(ZVOXDIM),...
' 1 32 32676 -32676 1 0 16'];
eval(str);
end

% **************
% Total Variance
fid=fopen(['TotalVariance.img'],'w',ENDIANOUT);
fwrite(fid,T,'float32');
fclose(fid);
str=['!makeheader4 TotalVariance.hdr ', num2str(XDIM), ' ', num2str(YDIM), ' ', num2str(ZDIM),...
' ', num2str(XVOXDIM), ' ', num2str(YVOXDIM), ' ', num2str(ZVOXDIM),...
' 1 32 32676 0 1 0 16'];
eval(str);

% ****************************
% Reduced Model Total Variance
fid=fopen(['ReducedModelVariance.img'],'w',ENDIANOUT);
fwrite(fid,EplusH_reduced,'float32');
fclose(fid);
str=['!makeheader4 ReducedModelVariance.hdr ', num2str(XDIM), ' ', num2str(YDIM), ' ', num2str(ZDIM),...
' ', num2str(XVOXDIM), ' ', num2str(YVOXDIM), ' ', num2str(ZVOXDIM),...
' 1 32 32676 0 1 0 16'];
eval(str);

% *************
% Wilk's Lambda
fid=fopen(['Lambda.img'],'w',ENDIANOUT);
fwrite(fid,Lambda,'float32');
fclose(fid);
str=['!makeheader4 Lambda.hdr ', num2str(XDIM), ' ', num2str(YDIM), ' ', num2str(ZDIM),...
' ', num2str(XVOXDIM), ' ', num2str(YVOXDIM), ' ', num2str(ZVOXDIM),...
' 1 32 1.0 0.0 1 0 16'];
eval(str);

% *********
% Omnibus F
fid=fopen(['F_' num2str(VARIABLES) '_' num2str(N-VARIABLES-1) '.img'],'w',ENDIANOUT);
fwrite(fid,F_omnibus,'float32');
fclose(fid);
str=['!makeheader4 ', ['F_' num2str(VARIABLES) '_' num2str(N-VARIABLES-1) '.hdr '],...
num2str(XDIM), ' ', num2str(YDIM), ' ', num2str(ZDIM),...
' ', num2str(XVOXDIM), ' ', num2str(YVOXDIM), ' ', num2str(ZVOXDIM),...
' 1 32 32676 0 1 0 16'];
eval(str);

% ******************************************************
% Wilks Reduced Model: Corrected for Confounding Effects
fid=fopen(['Lambda_reduced.img'],'w',ENDIANOUT);
fwrite(fid,Lambda_reduced,'float32');
fclose(fid);
str=['!makeheader4 ', ['Lambda_reduced.hdr '],...
num2str(XDIM), ' ', num2str(YDIM), ' ', num2str(ZDIM),...
' ', num2str(XVOXDIM), ' ', num2str(YVOXDIM), ' ', num2str(ZVOXDIM),...
' 1 32 1.0 0.0 1 0 16'];
eval(str);

% **************************************************
% F Reduced Model: Corrected for Confounding Effects
fid=fopen(['F_reduced_' num2str(VARIABLES-sum(S)) '_' num2str(N-VARIABLES-1) '.img'],'w',ENDIANOUT);
fwrite(fid,F_reduced,'float32');
fclose(fid);
str=['!makeheader4 ', ['F_reduced_' num2str(VARIABLES-sum(S)) '_' num2str(N-VARIABLES-1) '.hdr '],...
num2str(XDIM), ' ', num2str(YDIM), ' ', num2str(ZDIM),...
' ', num2str(XVOXDIM), ' ', num2str(YVOXDIM), ' ', num2str(ZVOXDIM),...
' 1 32 32676 0 1 0 16'];
eval(str);

% ***************************************
% F Reduced Model: Thresholded Mask Image
fid=fopen(['F_reduced_mask.img'],'w',ENDIANOUT);
  fwrite(fid,F_reduced_mask,'int8');
  fclose(fid);
  str=['!makeheader4 F_reduced_mask.hdr ', num2str(XDIM), ' ', num2str(YDIM), ' ', num2str(ZDIM),...
  ' ', num2str(XVOXDIM), ' ', num2str(YVOXDIM), ' ', num2str(ZVOXDIM), ' 1 8 1 0 1 0 2'];
  eval(str);

% *******************************************
% Individual Regressor and CONTRASTS F Values
for j=1:CONTRASTS,
   fid=fopen(['F_regressor_' num2str(j) '_' num2str(1) '_' num2str(N-VARIABLES-1) '.img'],'w',ENDIANOUT);
   fwrite(fid,F_indiv(j,:),'float32');
   fclose(fid);
   str=['!makeheader4 ', ['F_regressor_' num2str(j) '_' num2str(1) '_' num2str(N-VARIABLES-1) '.hdr '],...
num2str(XDIM), ' ', num2str(YDIM), ' ', num2str(ZDIM),...
' ', num2str(XVOXDIM), ' ', num2str(YVOXDIM), ' ', num2str(ZVOXDIM),...
' 1 32 32676 0 1 0 16'];
   eval(str);
end

% ***************************************************
% For Display, Plot Results of Maximum Omnibus F-test
MaxF_reduced=find(F_reduced==max(F_reduced));
z=floor(MaxF_reduced./(XDIM*YDIM))+1;
x(:,:)=Z(:,:,z);
yprime=x'*Beta(:,MaxF_reduced);

figure;
y=Y(:,MaxF_reduced)-mean(Y(:,MaxF_reduced));
plot(y,'r.');
hold
plot(yprime,'b');
hold off
xlabel('Time');
ylabel('Normalized MR Units');
title(['Overall Fitted Model at Voxel ' num2str(MaxF_reduced) ' on slice ' num2str(z)]);

% **********************************
% Generate Plot of Corrected Results
Ycorrected=(Y(:,MaxF_reduced)-mean(Y(:,MaxF_reduced)))-x'*(Beta(:,MaxF_reduced).*S');
yprime=x'*Beta(:,MaxF_reduced)-x'*(Beta(:,MaxF_reduced).*S');

figure;
plot(Ycorrected,'r.');
hold
plot(yprime,'b');
hold off
title(['MR Signal Corrected for Run to Run Effects at Voxel ' num2str(MaxF_reduced) ' on slice ' num2str(z)]);
xlabel('Time');
ylabel('Normalized MR Units');


% *******************************************************************************
% Find Local Maxima for Reduced Model Image and Write Results Tables to Text File
% and Time Series to a Matlab file

if(LOCAL_MAX_PLOTS==1)
 !rm LocalMaxPlots.eps;
end

outfile='TrackingLocalMax.txt';
eval(['!rm ' outfile]);
disp(' ');
disp('Now computing parameter image FWHM and identifiying local maxima...');

G=[];
R=[];
prob=1.0;

for i=1:length(S),
   if(S(i)==1)
      [Sx, Sy, Sz, U, Resels]=ComputeResels2(Beta(i,:), mask, XDIM, YDIM, ZDIM);
      G=[G;[Sx,Sy,Sz]];
      R=[R;Resels];
   end
end
Sxmean=mean(G(:,1));
Symean=mean(G(:,2));
Szmean=mean(G(:,3));

R=length(V)./((Sxmean*Symean*Szmean).^(2/3));
%B=(4/3)*pi*sqrt(det(S./2));
Reselsmean=mean(R);
LocalMax=FindLocalMax(F_reduced, XDIM, YDIM, ZDIM, Sxmean, Symean, Szmean, Reselsmean, Fcrit);
disp([num2str(length(LocalMax)) ' local maxima identified']);
PeakReducedModelTimeseries=zeros(N,length(LocalMax));

if(LOCAL_MAX_PLOTS==1)  % Make a new figure if local max plots are desired
         figure;
         !rm LocalMaxPlots.*
end

[P Q]=size(LocalMax);

for i=1:P,

      xlocation=LocalMax(i,2);
      ylocation=LocalMax(i,3);
      Index=LocalMax(i,5);
      %zlocation=floor(Index./(XDIM*YDIM))+1;
      zlocation=LocalMax(i,4);
      
      ANOVATables(outfile,xlocation,ylocation,zlocation,...
      N,VARIABLES,Beta(:,Index),T(Index),...
      F_omnibus(Index),(VARIABLES-sum(S)),F_reduced(Index),...
      Lambda(Index),Lambda_reduced(Index),EplusH_reduced,F_indiv(:,Index));
      
      x=Z(:,:,zlocation);
      yprime=x'*(Beta(:,Index).*S');
      y=Y(:,Index)-mean(Y(:,Index));
      PeakReducedModelTimeseries(:,i)=(y-yprime);
      yint=x'*(Beta(:,Index).*(1-S)');

      if(LOCAL_MAX_PLOTS==1)
         
         subplot(1,2,1);
         plot((y-yprime),'r.');
         hold on;
         plot(yint,'b');
         xlabel('Time');
         ylabel('Y');
         title(['Local Max F Reduced ' num2str(i) ' = ' num2str(F_reduced(Index)) ' (Lambda Reduced = ' num2str(Lambda_reduced(Index)) ') on voxel ' num2str([xlocation ylocation zlocation])]);
         hold off;
         %drawnow;
         
         subplot(1,2,2);
         plot((y-yprime),yint,'g.');
         %hold on;
         %b=linspace(min(y-yprime),max(y-yprime),10);
         %plot(b,b,'r');
	 xlabel('Y');
         ylabel('Y Predicted');
         title(['R-sqr = ' num2str(1-Lambda_reduced(Index))]);
         axis square;
         hold off;
         drawnow;
         %pause(0.5);
         print LocalMaxPlots.eps -dpsc -append;
      end

end

B=Beta(:,LocalMax(:,5));
save PeakReducedModelTimeseries PeakReducedModelTimeseries LocalMax Z B C;

% ***************************************************************************
% Fill Structure for Christian Buechel's Structural Equation Modeling Routine
%
% Struct array Data
% -----------------
% Data(k).X      - n x m matrix of time-series of m regions
% Data(k).L      - 3 x m matrix for centers of VOI for m regions
% Data(k).Rname	 - m x l matrix containing string descriptors (length l) 
%		   for m regions 
% Data(k).useit	 - index vector, indicating which regions of X are used 
%		   (max(useit) <= m)
% Data(k).FlagOb - Vector of size(useit) indicating whether a variable 
%		   is observed (=1) or latent (=0)
% Data(k).mod	 - Flag to indicate whether backwards modulatory connections 
%		   should be introduced for all (driving) connections
% Data(k).RT	 - Repetition time for experiment (to correct for 
%		   autocorrelation)
% Data(k).U 	 - 2 x g matrix coding g unilateral connections 
%		   (U(1,:) = sources and U(2,:) = destination)
% Data(k).B 	 - 2 x f matrix coding f bidirectional connections 
%		   between (B(1,:) and B(2,:)
%
%
% Cell array C
% --------------
% C{h}(f)	 - struct array cantaining f constraints
% C{h}(f).value	 - constrain the following paths to value, 
%		   if value = Inf impose equality constraint
% C{h}(f).conn	 - 2 x n matrix, (1,n) is index into 'Data(k)', 
%		   (2,n) denotes which paths within "Data"
%		   !note, counting starts with U and then B	 
%
% Struct Misc
% -----------
% Misc.Output 	 - Degree of output wanted: 0 = none, 1 = some, 2 = all tables, 		   
% Misc.Descr     - String to describe the analysis
% Misc.random	 - Use random starting estimates
%
%
% OUTPUT:
%
% Struct array G
% --------------
% G(h).value	 - constrain the following paths to value, 
%		   if value = Inf impose equality constraint
% G(h).conn	 - 2 x n matrix, (1,n) is index into Data array (k),
%		   (2,n) denotes which paths within "Data"
%		   !note, counting starts with U and then B	 
% G(h).Free      - Matrix of free parameters
% G(h).xfix      - fixed parameters
% G(h).SEM	 - see SEM
% G(h).x0        - starting parameters for optimisation 
% G(h).chi_sq	 - chi squared goodness of fit
% G(h).df	 - degrees of freedom
% G(h).p	 - p value
% G(h).RMSEA	 - RMSEA fit index
% G(h).ECVI	 - scale AIC fit index 
%
% Struct array SEM 
% ----------------
% SEM(k).Ustring	- Strings describing uni-directional connections
% SEM(k).Bstring	- Strings describing bi-directional connections
% SEM(k).Rname		- Strings describing regions (only updated 
%		          if modulatory connections are modelled)
% SEM(k).Constring	- Strings describing all
% SEM(k).Fil	 	- Filter to calculate covariances with latent variables 
% SEM(k).Cov		- covariance matrix
% SEM(k).df     	- degrees of freedom
% SEM(k).ConX	 	- matrix of unidirectional paths
% SEM(k).ConZ	 	- matrix of bidirectional paths
% SEM(k).A		- asymmetric path coefficients
% SEM(k).S		- symmetric path coefficients
% SEM(k).AL		- Modification indices for asymmetric path coefficients (*)
% SEM(k).SL		- Modification indices for symmetric path coefficients
% SEM(k).f		- fit index from optimisation
% SEM(k).Res		- residual covariances
% SEM(k).Est 		- estimated covariances
% SEM(k).L 		- modification indices for each parameter of
%
if (SEM==1)

Rname=[];
for i=1:P,
  if(i<10)
   reg=['00' num2str(i)];
  elseif((i>9)&(i<100))
   reg=['0' num2str(i)];
  elseif((i>99)&(i<1000))
   reg=[num2str(i)];
  end
  Rname=[Rname; reg];
end

     st = struct('X', PeakReducedModelTimeseries(:,1:P)',...
		 'L', LocalMax(:,2:4),...
		 'Rname',Rname,...
		 'useit', [1:P],...
		 'FlagOb', ones(size([1:P])),...
		 'mod', 1,...
		 'RT', TR,...
                 'U', [1 2;...
		      1 2],...         	
		 'B', [1 2;...
	               1 2]);
     C1 =   struct('value',{[]},...
		  'conn',{[]});
     C2 =   struct(    'value',{Inf,-.5},...
		  'conn',{[1 2;...
		  	   1 1],...
			   [3;...
			    1]});
     clear C;
     Data(1)=st;
     C{1}=C1;
     C{2}=C2;
     Misc =   struct('Output',2,...
		  'Descr' ,'Reduced Model Regression Local Max Data',...
		  'random',0);

     save SEM Data C Misc;

end

% *******
% The End
disp('All Done.');
disp(['Total Program Running Time = ' num2str(floor(toc/60)) ':' num2str(60*mod(toc/60,1)) ]);




% %%%%%%%%% %
% JVH, 2001 %
% %%%%%%%%% %







