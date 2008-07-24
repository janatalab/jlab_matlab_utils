function [Sx, Sy, Sz, U, Resels]=ComputeResels2(X, Mask, XDIM, YDIM, ZDIM)
% Compute FWHM of the 3D Autocorrelation Function
% and use to Estimate the Number of Resolvable Elements
% in the Image.

DIMS=[XDIM YDIM ZDIM];

%fid=fopen(ImageName,'r');
%  X=fread(fid,'float32');
%fclose(fid);

%fid=fopen('Mask.img','r');
%  Mask=fread(fid);
%fclose(fid);
M=sum(Mask);

Z=reshape(X,XDIM,YDIM,ZDIM);

% Determine 3D Autocorrelation Function of the Image
z=fftn(Z);
w=z.*conj(z);
W=real(ifftn(w));
W=W./W(1,1,1);  % So that the max is unity
W=fftshift(W);  % So that the DC offset is a the center of the volume

S=zeros(3,3);
SumW=0;

% Measure Upper Triangle Cross-Products of the 3D Autocorrelation Function
for x=1:DIMS(1),
	for y=1:DIMS(2),
		for z=1:DIMS(3),
			S(1,1)=S(1,1)+abs(W(x,y,z)).*(x-(DIMS(1)+1)./2).^2;
end
end
end

for x=1:DIMS(1),
	for y=1:DIMS(2),
		for z=1:DIMS(3),
			S(2,2)=S(2,2)+abs(W(x,y,z)).*(y-(DIMS(2)+1)./2).^2;
end
end
end

for x=1:DIMS(1),
	for y=1:DIMS(2),
		for z=1:DIMS(3),
			S(3,3)=S(3,3)+abs(W(x,y,z)).*(z-(DIMS(2)+1)./2).^2;
end
end
end

for x=1:DIMS(1),
	for y=1:DIMS(2),
		for z=1:DIMS(3),
			S(1,2)=S(1,2)+abs(W(x,y,z)).*(x-(DIMS(1)+1)./2).*(y-(DIMS(2)+1)./2);
end
end
end

for x=1:DIMS(1),
	for y=1:DIMS(2),
		for z=1:DIMS(3),
			S(1,3)=S(1,3)+abs(W(x,y,z)).*(x-(DIMS(1)+1)./2).*(z-(DIMS(3)+1)./2);
end
end
end

for x=1:DIMS(1),
	for y=1:DIMS(2),
		for z=1:DIMS(3),
			S(2,3)=S(2,3)+abs(W(x,y,z)).*(y-(DIMS(2)+1)./2).*(z-(DIMS(3)+1)./2);
end
end
end

% Compute Sum of Absolute Autocorrelation Values
SumW=sum(sum(sum(abs(W))));

% Reflect and Scale the Autocorrelation Cross-Product Matrix by 2*log(2)./SumW
%S=2*log(2)*(S+S'.*(1-eye(3,3)))./SumW;
S=(S+S'.*(1-eye(3,3)))./SumW;

% Obtain SVD Decomposition of S
[U V P]=svd(S);

% Estimated FWHM in Each Axis
Sx=sqrt(V(1,1));
Sy=sqrt(V(2,2));
Sz=sqrt(V(3,3));

% Compute the Area of an Ellipsoid With These Major Axes
A=(S(1,1)*S(2,2)*S(3,3))^(-1/6);
B=(4/3)*pi*sqrt(det(S./2));

% Contour Method: Find ISO contour at W==0.5 in X and Y
%C=contourc(W(:,:,(ZDIM+1)./2),[0.5 0.5]);
%F(1,:)=C(1,:)-(XDIM+1)./2;
%F(2,:)=C(2,:)-(YDIM+1)./2;
%b=F(:,2:end)*F(:,2:end)';
%Sx=sqrt(b(1,1));
%Sy=sqrt(b(2,2));

% Find ISO contour at W==0.5 in Z
%clear C F;
%a=permute(W,[3 1 2]);
%C=contourc(a(:,:,ceil((YDIM+1)./2)),[0.5 0.5]);
%F(1,:)=C(1,:)-(ZDIM+1)./2;
%b=F(:,2:end)*F(:,2:end)';
%Sz=sqrt(b(1,1));

% Compute Resolveable Elements 
%Resels=M.*A;
Resels=M.*B^(-1/3);

% Clear Temp Variables
clear X Z z w Mask M V a F;

return;







































