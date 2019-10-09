%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Distribution of proposed sampling points that 
%%% generate low coherence spherical harmonics sensing matrices
%%%
%%%  Created by Arya Bangun at TI RWTH Aachen 2019 09.10.2019



clear all
close all

%% Generate unit sphere
figure;
s1=subplot(1,3,1);
rIn=1;
phiRadOut=0:5:360;
thetaRadOut=0:5:180;
[phi,theta]=meshgrid(phiRadOut,thetaRadOut);
phi=deg2rad(phi);
theta=deg2rad(theta);
x=rIn*sin(theta).*cos(phi);
y=rIn*sin(theta).*sin(phi);
z=rIn*cos(theta);
mesh(x,y,z)
axis equal
grid on
box off
axis off
view([130,30])
colormap([0.5 .5 .5])
hold on

%% Load and plot sampling points
load SH_N1024
idx = [1 9 length(m)];
ang=total_angles.proposed{idx(1)};
theta_sampl = ang(:,1);
phi_sampl =  ang(:,2);
theta_sampl=repmat(theta_sampl,1,size(phi_sampl,2));
x_sampl=rIn*sin(theta_sampl).*cos(phi_sampl);
y_sampl=rIn*sin(theta_sampl).*sin(phi_sampl);
z_sampl=rIn*cos(theta_sampl);
X = [x_sampl(:)];
Y = [y_sampl(:)];
Z = [z_sampl(:)];
scatter3(X(:),Y(:),Z(:),10,'ko','MarkerFaceColor','r','LineWidth',2)
title(['m = ', num2str(m(idx(1)))]);

%% Proposed m=500 B=32

s2=subplot(1,3,2);
rIn=1;
phiRadOut=0:5:360;
thetaRadOut=0:5:180;
[phi,theta]=meshgrid(phiRadOut,thetaRadOut);
phi=deg2rad(phi);
theta=deg2rad(theta);
x=rIn*sin(theta).*cos(phi);
y=rIn*sin(theta).*sin(phi);
z=rIn*cos(theta);
mesh(x,y,z)
axis equal
grid on
box off
axis off
view([130,30])
colormap([0.5 .5 .5])
hold on
ang=total_angles.proposed{idx(2)};
theta_sampl = ang(:,1);
phi_sampl =  ang(:,2);
theta_sampl=repmat(theta_sampl,1,size(phi_sampl,2));
x_sampl=rIn*sin(theta_sampl).*cos(phi_sampl);
y_sampl=rIn*sin(theta_sampl).*sin(phi_sampl);
z_sampl=rIn*cos(theta_sampl);
X = [x_sampl(:)];
Y = [y_sampl(:)];
Z = [z_sampl(:)];
scatter3(X(:),Y(:),Z(:),10,'ko','MarkerFaceColor','r','LineWidth',2)
title(['m = ', num2str(m(idx(2)))]);

%% Proposed m=900 B=32

s3=subplot(1,3,3);
rIn=1;
phiRadOut=0:5:360;
thetaRadOut=0:5:180;
[phi,theta]=meshgrid(phiRadOut,thetaRadOut);
phi=deg2rad(phi);
theta=deg2rad(theta);
x=rIn*sin(theta).*cos(phi);
y=rIn*sin(theta).*sin(phi);
z=rIn*cos(theta);
mesh(x,y,z)
axis equal
grid on
box off
axis off
view([130,30])
colormap([0.5 .5 .5])
hold on
ang=total_angles.proposed{idx(end)};
theta_sampl = ang(:,1);
phi_sampl =  ang(:,2);
theta_sampl=repmat(theta_sampl,1,size(phi_sampl,2));
x_sampl=rIn*sin(theta_sampl).*cos(phi_sampl);
y_sampl=rIn*sin(theta_sampl).*sin(phi_sampl);
z_sampl=rIn*cos(theta_sampl);
X = [x_sampl(:)];
Y = [y_sampl(:)];
Z = [z_sampl(:)];
scatter3(X(:),Y(:),Z(:),10,'ko','MarkerFaceColor','r','LineWidth',2)
%set(gca, 'Position', [0.1 0.1 0.3 0.85])
title(['m = ', num2str(m(idx(3)))]);
