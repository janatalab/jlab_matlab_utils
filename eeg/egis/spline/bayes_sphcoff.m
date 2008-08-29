function [sp_coff,error_check, sp_mat] = bayes_sphcoff(nmax,xelec,yelec,zelec,v,lambda)
%function [sp_coff,error_check, sp_mat] = bayes_sphcoff(nmax,xelec,yelec,zelec,v,lambda)
%
%	calculates spherical harmonic coefficients based on jerry's notes on Bayesian formulation for linear systems
%
%	nmax = maximum order allowed for spherical harmonic sequences. 
%			recommended (maximum order for full rank model)
%						129 channel EGI system nmax = 8
%	xelec = x locations of electrode array
%	yelec = y locations of electrodes
%	zelec = z location of electrodes
%	v = data
%	lambda = smoothing parameter (one day this will become noise)
%
%	outputs
%	sp_coff = spherical harmonic coefficients
%	error_check = error estimate
%	sp_mat = spherical harmonic matrix
if ~(nargin == 5|nargin == 6)
	error('improper parameter list')
end;

[azelec,elelec,relec] = cart2sph(xelec,yelec,zelec);
elelec = (pi/2)*ones(size(elelec,1),size(elelec,2)) - elelec;

icount = 1;

for j=1:nmax
	sp_mat(icount:icount+2*j,:) = sph_elec(xelec,yelec,zelec,j);
	icount = icount+2*j+1;
end;

if icount~= (nmax+1).^2 
	error('Dope')
end;

ls_mat = sp_mat*sp_mat';

if nargin == 6
	for i = 1:size(ls_mat,1)
		ls_mat(i,i) = ls_mat(i,i) + lambda;
	end;
end;

lv_mat = sp_mat*v';

sp_coff = ls_mat\ lv_mat;

error = sp_mat'*sp_coff - v';

error_check(1) = sum(abs(error))/129;

