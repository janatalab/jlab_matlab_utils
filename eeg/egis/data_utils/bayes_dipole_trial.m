function [m,deviations] = bayes_dipole_trial(a,v,sigma_v,sigma_m)
%function [m,deviations] = bayes_dipole(a,v,sigma_v,sigma_m)
%
%   a = forward coefficients (from forward model)
%	v = data
%	sigma_v = std deviation of data
%   sigma_m = std deviation of moments (for prior)
%
%	outputs
%	m = dipole moments at maximum posterior likelihood
%	deviations = std deviations
%     deviations(1) = std deviation of potential errors 
%          (should be ~= sigma_v)
%     deviations(2) = std deviation of moments
%          (should be ~= sigma_m)

if ~(nargin == 4)
	error('improper parameter list')
end;


ls_mat = a'*a + ((sigma_v/sigma_m)^2  * eye(size(a,2)));

inv_ls_mat = inv(ls_mat);
lv_mat = a'*v;

m = inv_ls_mat*lv_mat;

error = a*m - v;

deviations(1,:) = sqrt(sum(error.^2)/size(error,1));

deviations(2,:) = sqrt(sum(m.^2)/size(m,1));
