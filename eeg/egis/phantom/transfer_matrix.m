function [potmat, lapmat, cortmat] = transfer_matrix(nmax,s12,s13,r1,r2,r3,r4,rz,sources,electrodes);
%[potmat,lapmat,cortmat] = transfer_matrix(nmax,s12,s13,r1,r2,r3,r4,rz,sources,electrodes);
%
[a,b,c,d,e,f,g] = fender_4(nmax,s12,s13,r1,r2,r3,r4,rz);
h = f+g;

for i = 1:nmax
l(i) = i*(i+1)*h(i)/r4/r4;
end;
for i = 1:size(a,2)
a(i) = a(i)+i*((rz/r1)^(i-1))/(4*pi*r1^2);
end;

cos_mat = ang_dist(sources,electrodes);

potmat = zeros(size(cos_mat));
lapmat = zeros(size(cos_mat));
cortmat = zeros(size(cos_mat));
for i = 1:size(cos_mat,1)
 	pmat = legendre_r(cos_mat(i,:),nmax);
	hmat = h'*ones(1,size(cos_mat,2));
	lmat = l'*ones(1,size(cos_mat,2));
	amat = a'*ones(1,size(cos_mat,2));
	hmat = hmat.*pmat;
	lmat = lmat.*pmat;
	amat = amat.*pmat;
	potmat(i,:) = sum(hmat);
	lapmat(i,:) = sum(lmat);
	cortmat(i,:) = sum(amat);
end;
potmat = potmat*4*pi*r4*r4;
lapmat = lapmat*4*pi*r4*r4;
cortmat = cortmat*4*pi*r4*r4;






	