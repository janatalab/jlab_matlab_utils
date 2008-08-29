function v_lap = spherical_lap(welec,x,y,z,xs,ys,zs,p,q)
%
% calculates laplacian on spherical surface
%
welec = welec.^2;
n = size(x,2);
ns = size(xs,2);
[azs,els,rs] = cart2sph(xs,ys,zs);
els = (pi/2)*ones(size(els,1),size(els,2)) - els;
[az,el,r] = cart2sph(x,y,z);
el = (pi/2)*ones(size(el,1),size(el,2)) - el;

% calculate trig functions to keep things easier

st = sin(els);
ct = cos(els);
sp = sin(azs);
cp = cos(azs);
s2t = sin(2*els);
c2t = cos(2*els);
s2p = sin(2*azs);
c2p = cos(2*azs);

% calculate utta uppa utt upp ucc

% the real enchilada is uuxyz 
%
% note for future:  this is where i have to rederive for arbitrary surface normals
%

uuxyz = 2*q(4)+2*q(6)+2*q(10) - (2*st.*(q(2)*cp+q(3)*sp)./rs +2*q(7)*ct./rs+6*(st.^2).*(q(4)*(cp.^2)+q(6)*(sp.^2)+q(5)*sp.*cp)+6*st.*ct.*(q(8)*cp+q(9)*sp)+6*q(10)*ct.^2);


%
% its the unfortunate time where we deal witht logarithmic mess
%

%
% zero registers
%

smpp = zeros(size(st,1),size(st,2));
smtt = zeros(size(st,1),size(st,2));
smp = zeros(size(st,1),size(st,2));
smt = zeros(size(st,1),size(st,2));
ttcomp = zeros(size(st,1),size(st,2));
rrcomp = zeros(size(st,1),size(st,2));

for j = 1:n
 	a = r(j)*(st.*cp-sin(el(j))*cos(az(j))*ones(size(st,1),size(st,2)));
	b = r(j)*(st.*sp-sin(el(j))*sin(az(j))*ones(size(st,1),size(st,2)));
	c = r(j)*(ct-cos(el(j))*ones(size(st,1),size(st,2)));
	str = a.^2+b.^2+c.^2;
	strw = str+welec*ones(size(st,1),size(st,2));
	comterm = 4*str./strw-((str./strw).^2)+2*log(strw);
	comterm2 = 2*(2*str.*log(strw)+(str.^2)./strw);
	tcomp = 3*comterm2+4*str.*comterm;
	dr = 2*(a.*st.*cp+b.*st.*sp+c.*ct);
	d2r2 = 2;
	rcomp = dr.*comterm2+d2r2*r(j)*comterm2/2+r(j)*(dr.^2).*comterm;
	ttcomp = ttcomp + p(j)*tcomp;
	rrcomp = rrcomp + p(j)*rcomp/r(j);
end;
v_lap = - (ttcomp+uuxyz-rrcomp);










