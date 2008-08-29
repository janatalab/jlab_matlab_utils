function [interpdata, Hpatch] = topo3D(data,badchan,minmax,sigvals,thetitle);
%topo3D(data,badchan,minmax,sigvals);
%makes a 3D topographic map.  Plots 2 standard views in separate figure
%windows.  Each of these views can be rotated by the controls in the figure
%window.  
%data = 1by129 data at the EEG electrode positions
%badchan = list of bad channels
%minmax = [minimum maximum] data range for plots.  If an empty matrix 
%is provided the program automatically sets the range
%sigvals = (optional) 1by129 array of p values. if p < 0.05 the color of the electrode
%positions is set to green otherwise its gray.  This parameter can be skipped.  
global hotncold hotter
global stimchan EYE
figure(100)
hotncold2 = hotncold;
hotter2 = jet;
hotncold2(1,:) = 0.5;
hotter2(1,:) = 0.5;
close(100)
if nargin ~=3  & nargin ~=4 & nargin ~= 5 
    error('incorrect number of inputs');
end;
global ENVIRONMENT
load headmodel.mat
chans = [1:129];
chans = removebadchan(chans,badchan);
if nargin == 3
    sigvals = 0.01*ones(1,129);
    thetitle = [];
    if ~strcmp(ENVIRONMENT,'Lausanne')
    sigvals([stimchan EYE]) = 1;    
    else
    sigvals(EYE) = 1;
    end;
end;
if nargin == 4
thetitle = [];
end;
if isempty(minmax);
    if min(data(chans)) >= 0
        maxdata = max(data(chans));
        mindata = 0;
    else
        maxdata = max(abs(data(chans)));
        mindata = - maxdata;
    end
else
    mindata = minmax(1);
    maxdata = minmax(2);
end;

x = elp.x(chans)';
y = elp.y(chans)';
z = elp.z(chans)';
v = data(chans);
if mindata == 0
    v = round((v/maxdata+1/63)*63);
else
    v = round((v/maxdata+32/63)*63);
end;
list1 = find(vertex_matrix(:,3) > 0.02 & vertex_matrix(:,1) > -0.04);
list2 = find(vertex_matrix(:,3) > 0.00 & vertex_matrix(:,1) <= -0.04);
xs = vertex_matrix([list1; list2],1);
ys = vertex_matrix([list1; list2],2);
zs = vertex_matrix([list1; list2],3);
w= 10^-15;
[k,kinv,a,ainv,e]= k_and_e(w,x,y,z);
[p,q,error_check]= mateqs(w,x,y,z,v,k,kinv,a,ainv,e);
interp3D = interp_3d(w,x,y,z,xs,ys,zs,p,q);
interpdata = ones(length(vertex_matrix(:,1)),1); 
interp3D(find(interp3D <= 2)) = 2;
interp3D(find(interp3D >= 64)) = 64;
interpdata(list1) = interp3D(1:length(list1));
interpdata(list2) = interp3D(length(list1)+1:end);
figure
Hpatch = patch('Vertices',vertex_matrix,'Faces',face_matrix,...
'FaceVertexCData',interpdata,'facecolor','interp','edgecolor','k','CDataMapping','direct');
if mindata < 0 
    colormap(hotncold2)
else
    colormap(hotter2)
end
hold on 
global stimchan EYE
chanselec = [1:129];
chanselec = removebadchan(chanselec,[stimchan EYE]);
h = plot3(elp.x(chanselec),elp.y(chanselec),elp.z(chanselec),'ko');
set(h,'MarkerFaceColor',[0.5 0.5 0.5]);
set(h,'MarkerSize',4);
sigchans = find(sigvals < 0.05);
sigchans = removebadchan(sigchans, badchan);
h = plot3(elp.x(sigchans),elp.y(sigchans),elp.z(sigchans),'ko');
set(h,'MarkerFaceColor','k');
set(h,'MarkerSize',9);
axis image
axis off
h = colorbar;
set(h,'YTick',[1 33 65]);
set(h,'YTickLabel',[minmax(1) sum(minmax)/2 minmax(2)]');
axis off
view(-90,50)
title(thetitle);
%
