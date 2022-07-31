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
thetaInDeg = (-90:deltaTheta:90)';
thetaIn = thetaInDeg/180*pi;

thetaOutDeg = (-90:deltaTheta:90)';
thetaOut = thetaOutDeg/180*pi;

[thetaInDegMesh, thetaOutDegMesh] = meshgrid(thetaInDeg, thetaOutDeg);
[thetaInMesh, thetaOutMesh] = meshgrid(thetaIn, thetaOut);

thetaIn = reshape(thetaInMesh, [], 1);
thetaOut = reshape(thetaOutMesh, [], 1);

%%
Sa = sinc(b./lambda*(sin(thetaOut)+sin(thetaIn)));

RCS = numUnits*4*pi*(abs(C))^2*(a*b)^2.*(cos(thetaIn)).^2.*Sa;

RCS = reshape(RCS,size(thetaInMesh));

%%
figure
mesh(thetaInDegMesh, thetaOutDegMesh, RCS)

xlabel('$\theta^i$ (degree)','interpreter','latex')
ylabel('$\theta^s$ (degree)','interpreter','latex')
zlabel('$E [ \sigma(\theta^s; \theta^i) ]$','interpreter','latex')
axis tight
rotate3d on
% axis equal
% view(45, 30);
% view(0, 90);
box on

xh = get(gca,'XLabel'); % Handle of the x label
set(xh, 'Units', 'Normalized')
pos = get(xh, 'Position');
set(xh, 'Position',pos.*[1.1,-0.5,1],'Rotation',15)
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[0.7,-0.6,1],'Rotation',-25)

exportgraphics(gcf, 'ArrayRCS_RandomPhaseShift.pdf');