clear all;close all;clc;
%%
deltaTheta = 0.1;

C = -2;

lambda = 1;

%%
a = 0.4;
b = 0.1;

d = 0.45;

numUnits = 30;
unitsPos = (0:1:numUnits-1)*d;

%%
anglesIncident = [10, 70]/180*pi;
amplitudeIncident = [1, 2];

anglesSteer = -50/180*pi;

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
w = ones(numUnits,1);
wTemp = exp(-1i*2*pi.*sin(anglesIncident(1)).*unitsPos).*exp(-1i*2*pi.*sin(anglesSteer)*unitsPos);
w = wTemp.';

Response = C*a*b.*exp(1i*2*pi.*(kron(unitsPos, sin(thetaIn))+kron(unitsPos, sin(thetaOut))))*w;

Response = reshape(abs(Response),size(thetaInMesh));

%%
figure
mesh(thetaInDegMesh, thetaOutDegMesh, abs(Response).^2)

xlabel('$\theta^i$ (degree)','interpreter','latex')
ylabel('$\theta^s$ (degree)','interpreter','latex')
zlabel('$\sigma(\theta^s; \theta^i)$','interpreter','latex')
axis tight
rotate3d on
% axis equal
view(45, 30);
% view(0, 90);
box on