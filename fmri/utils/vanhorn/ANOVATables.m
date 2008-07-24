function ANOVATables(outfile, x, y, z, N, p, Beta, TotalVariance, Fomni, p_int, FReduced, Lambda, LambdaReduced, EplusHReduced,F_contrast,C,C_names)
% ANOVATables(outfile, x, y, z, N, p, Beta, TotalVariance, Fomni, p_int,
%             FReduced, Lambda, LambdaReduced, EplusHReduced,F_regressor,C,C_names) 

% 03/24/02 PJ  Fixed handling of contrasts involving more than one regressor.

% Original script written by Jack Van Horn.

mode='a';

fid=fopen(outfile,mode);

% Omnibus ANOVA Table
ndf=p;
ddf=(N-p-1);
prob=betainc((ddf./(ddf+ndf.*Fomni)),ddf/2,ndf/2);
fprintf(fid,'*************************************************************************\n');
fprintf(fid,'*** Voxel %i, %i, %i Omnibus Analysis of Variance for Regression Table ***\n',x,y,z);
fprintf(fid,'*************************************************************************\n\n');
fprintf(fid,'-------------------------------------------------------------------------\n');
fprintf(fid,'Source\t\tdf\t\tSS\t\tF\t\tp\n');
fprintf(fid,'-------------------------------------------------------------------------\n');
fprintf(fid,'Regression\t%i\t\t%8.4f\t%8.4f\t\t%6.5f\n', p, (1-Lambda)*TotalVariance, Fomni, prob);
fprintf(fid,'Error\t\t%i\t\t%8.4f\n', (N-p-1), Lambda*TotalVariance);
fprintf(fid,'-------------------------------------------------------------------------\n');
fprintf(fid,'Total\t\t%i\t\t%8.4f\n', (N-1), TotalVariance);
fprintf(fid,'-------------------------------------------------------------------------\n');
fprintf(fid,'R-Squared = %f\n',1-Lambda);
fprintf(fid,'R-Squared Adjusted for Shrinkage= %f\n\n',1-Lambda*(N-1)/(N-p-1));

% Reduced Model ANOVA Table
ndf=p_int;
ddf=(N-p-1);
prob=betainc((ddf./(ddf+ndf.*FReduced)),ddf/2,ndf/2);
fprintf(fid,'*******************************************************************************\n');
fprintf(fid,'*** Voxel %i, %i, %i Reduced Model Analysis of Variance for Regression Table ***\n',x,y,z);
fprintf(fid,'*******************************************************************************\n\n');
fprintf(fid,'-------------------------------------------------------------------------\n');
fprintf(fid,'Source\t\tdf\t\tSS\t\tF\t\tp\n');
fprintf(fid,'-------------------------------------------------------------------------\n');
fprintf(fid,'Regression\t%i\t\t%8.4f\t%8.4f\t\t%6.5f\n', p_int, (1-LambdaReduced)*EplusHReduced, FReduced, prob);
fprintf(fid,'Error\t\t%i\t\t%8.4f\n', (N-p-1), LambdaReduced*EplusHReduced);
fprintf(fid,'-------------------------------------------------------------------------\n');
fprintf(fid,'Total\t\t%i\t\t%8.4f\n', (N-1), EplusHReduced);
fprintf(fid,'-------------------------------------------------------------------------\n');
fprintf(fid,'Reduced Model R-Squared = %f\n', 1-LambdaReduced);
fprintf(fid,'Reduced Model R-Squared Adjusted for Shrinkage= %f\n\n', 1-LambdaReduced*(N-1)/(N-p-1));
fprintf(fid,'-------------------------------------------------------------------------\n');
fprintf(fid,'Regression Parameters for Variables Remaining in Model:\n');

% Individual contrasts
ncont = size(C,1);
if ~exist('C_names')
  [C_names(1:ncont)] = deal({''});
end

for icont=1:ncont
  ndf=p-sum(C(icont,:));
  ddf=(N-p-1);
  prob=betainc((ddf./(ddf+ndf.*F_contrast(icont))),ddf/2,ndf/2);
  
  % If this is a test of a single regressor, show the beta value
  if ndf == 1
    reg_idx = find(~C(icont,:));
    fprintf(fid,'Model: %s\tRegressor %i:\tBeta = %f, F(1,%i) = %8.4f, p = %f\n', C_names{icont}, reg_idx, Beta(reg_idx), (N-p-1), F_contrast(icont), prob);
  else
    fprintf(fid,'Model: %s\t, F(%d,%d) = %8.4f, p = %f\n', C_names{icont}, ndf, ddf, F_contrast(icont), prob);
  end
end
fprintf(fid,'\n\n#########################################################################\n\n\n');

fclose(fid);
