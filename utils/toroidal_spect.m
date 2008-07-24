function [amp, theta, phi, names, Rsqr, X] = toroidal_spect(data, dims,m, n)
% [amp, theta, phi] = toroidal_spect(data, dims);
%
% Computes coefficients using least-squares of terms for spherical harmonics of
% a function on a toroidal surface. 
%
% f(theta,phi) = acc(m,n)*cos(m*theta)*cos(n*phi) +
%                acs(m,n)*cos(m*theta)*sin(n*phi) +
%                asc(m,n)*sin(m*theta)*cos(n*phi) +
%                ass(m,n)*sin(m*theta)*sin(n*phi) 
%         
%  NOTE: Each term is a summation over 1:m,1:n, where m and n are the order of
%        the highest spherical harmonic along each dimension.
%
%  data: the surface function, arranged as a row vector with dimensions of
%  theta and phi given in dims. If data is a matrix, rows are individual
%  surfaces that are to be modelled.
%

% 03/25/02 Petr Janata

nunits_theta = dims(1);
nunits_phi = dims(2);
nunits_total = nunits_theta * nunits_phi;

if any(size(data) ==  1)
  data = data(:)';
  nsamps = 1;
  if length(data) ~= nunits_total
    error('toroidal_spect: length of data does not match product of dimensions')
  end
else
  nsamps = size(data,1);
  if size(data,2) ~= nunits_total
    error('toroidal_spect: data size dimension mismatch')
  end
end


%
% Compute each of the sin and cos components of the expansion
%
clear cos_theta sin_theta cos_phi sin_phi

theta = 2*pi*(0:nunits_theta-1)/nunits_theta;
phi = 2*pi*(0:nunits_phi-1)/nunits_phi;

for iharm = 1:max([m n])
  harm = iharm-1;
  tmp = ...
      repmat(cos(harm*theta)',1,nunits_phi);
  cos_theta(:,iharm) = tmp(:);
  
  tmp = ...
      repmat(sin(harm*theta)',1,nunits_phi);
  sin_theta(:,iharm) = tmp(:);
  
  tmp = repmat(cos(harm*phi),nunits_theta,1);
  cos_phi(:,iharm) = tmp(:);
  
  tmp = repmat(sin(harm*phi),nunits_theta,1);
  sin_phi(:,iharm) = tmp(:);
end

%
% Compute the terms of the expansion for which amplitude parameters will be estimated
%
clear coscos cossin sincos sinsin
for mharm = 1:m
  for nharm = 1:n
    coscos(mharm,nharm,:) = cos_theta(:,mharm) .* cos_phi(:,nharm);
    coscos_names{mharm,nharm} = sprintf('cc%d%d',mharm-1,nharm-1);
      
    cossin(mharm,nharm,:) = cos_theta(:,mharm) .* sin_phi(:,nharm);
    cossin_names{mharm,nharm} = sprintf('cs%d%d',mharm-1,nharm-1);

    sincos(mharm,nharm,:) = sin_theta(:,mharm) .* cos_phi(:,nharm);
    sincos_names{mharm,nharm} = sprintf('sc%d%d',mharm-1,nharm-1);

    sinsin(mharm,nharm,:) = sin_theta(:,mharm) .* sin_phi(:,nharm);
    sinsin_names{mharm,nharm} = sprintf('ss%d%d',mharm-1,nharm-1);
  end % for jharm=
end % for iharm=

% Reshape the matrices so it is easy to access the whole set of harmonics for
% any given output unit on the map.
coscos = reshape(coscos,m*n,nunits_total)';
cossin = reshape(cossin,m*n,nunits_total)';
sincos = reshape(sincos,m*n,nunits_total)';
sinsin = reshape(sinsin,m*n,nunits_total)';

% Save names containing indices of each of the terms
names = [coscos_names(:)' cossin_names(:)' sincos_names(:)' sinsin_names(:)'];

% Put them together into a single regression matrix
X = [coscos cossin sincos sinsin];

% Compute some costly matrices right now
pinvX = pinv(X'*X);
pinvXX = pinvX*X';

amp = zeros(nsamps,m*n*4);

% Estimate regression coefficients for each of the samples
for isamp = 1:nsamps
  y = data(isamp,:)';
  y = y-mean(y);
  b = pinvXX*y;
  y_prime = X*b;
  y_err = y - y_prime;
  Rsqr(isamp) = 1- (y_err'*y_err)/(y'*y);
  amp(isamp,:) = b';
end