clear all;close all;clc;
%%
deltaTheta = 0.1;

C = -2*1i;

lambda = 1;

%%
a = 0.1;
b = 0.1;

d = 0.5;

numUnits = 100;
unitsPos = (0:1:numUnits-1)*d;

%%
anglesIncident = 30/180*pi;

anglesSteer = -50/180*pi;
% anglesSteer = -30/180*pi;


%%
thetaInDeg = (-90:deltaTheta:90)';
thetaIn = thetaInDeg/180*pi;

thetaOutDeg = (-90:deltaTheta:90)';
thetaOut = thetaOutDeg/180*pi;

[thetaInDegMesh, thetaOutDegMesh] = meshgrid(thetaInDeg, thetaOutDeg);
[thetaInMesh, thetaOutMesh] = meshgrid(thetaIn, thetaOut);

thetaIn = reshape(thetaInMesh, [], 1);
thetaOut = reshape(thetaOutMesh, [], 1);

%%
wTemp = exp(-1i*2*pi.*sin(anglesIncident(1)).*unitsPos).*exp(-1i*2*pi.*sin(anglesSteer)*unitsPos);
w = wTemp.';

Transfer = abs(exp(1i*2*pi.*(kron(unitsPos, sin(thetaIn))+kron(unitsPos, sin(thetaOut))))*w);

RCS = 4*pi*(abs(C))^2*(a*b)^2.*(cos(thetaIn)).^2.*Transfer.^2;


Transfer = reshape(Transfer,size(thetaInMesh));

RCS = reshape(RCS,size(thetaInMesh));



%%
figure
mesh(thetaInDegMesh, thetaOutDegMesh, Transfer)

xlabel('$\theta^i$ (degree)','interpreter','latex')
ylabel('$\theta^s$ (degree)','interpreter','latex')
zlabel('$T(\theta^s; \theta^i)$','interpreter','latex')
axis tight
% rotate3d on
% axis equal
view(35, 50);
% view(0, 90);
box on

xh = get(gca,'XLabel'); % Handle of the x label
set(xh, 'Units', 'Normalized')
pos = get(xh, 'Position');
set(xh, 'Position',pos.*[1,1,1],'Rotation',-12)
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[1,1,1],'Rotation',34)

exportgraphics(gcf, 'Transfer_PhaseCompensation.pdf');


%%
figure
mesh(thetaInDegMesh, thetaOutDegMesh, RCS)

xlabel('$\theta^i$ (degree)','interpreter','latex')
ylabel('$\theta^s$ (degree)','interpreter','latex')
zlabel('$\sigma(\theta^s; \theta^i)$','interpreter','latex')
axis tight
% rotate3d on
% axis equal
view(35, 50);
% view(0, 90);
box on

xh = get(gca,'XLabel'); % Handle of the x label
set(xh, 'Units', 'Normalized')
pos = get(xh, 'Position');
set(xh, 'Position',pos.*[1,1,1],'Rotation',-12)
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[1,1,1],'Rotation',34)

exportgraphics(gcf, 'RCS_PhaseCompensation.pdf');