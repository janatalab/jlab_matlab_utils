function [dist,bcohpotential,bcohlaplacian] = baseline_coherence(s12,s13,r1,r2,r3,r4,rz);
[a,b,c,d,e,f,g] = fender_4(500,s12,s13,r1,r2,r3,r4,rz);
ang = ([0:33]/33)*27/9.2;
dist = ang*9.2;
z = cos(ang);
p = legendre_r(z,500);
h = (f+g)';
numpot = p.*((h.^2)*ones(1,size(p,2)));
numpotsum = sum(numpot);
denpot = ((h.^2)*ones(1,size(p,2)));
denpotsum = sum(denpot);
bcohpotential = (numpotsum./denpotsum).^2;
numlap =   p.*((h.^2.*[1:500]'.*[2:501]')*ones(1,size(p,2)));
numlapsum = sum(numlap);
denlap = ((h.^2.*[1:500]'.*[2:501]')*ones(1,size(p,2)));
denlapsum = sum(denlap);
bcohlaplacian =  (numlapsum./denlapsum).^2;





