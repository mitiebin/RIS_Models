clear all;close all;clc;

%%
speedLight = 299792458;

f = 3.3e9;
% f = 0.3e9;

lambda = speedLight/f;

% a = 0.45*lambda;
% b = 0.45*lambda;

a = 0.042;
b = 0.042;

thetaIncidentDeg = 0;
phiIncidentDeg = 0;

Gamma = -1;

%%
deltaTheta = 0.2;
deltaPhi = 0.2;

csvFilePath = 'realCell_3.3G_0.2Deg.csv';
% csvFilePath = 'rcsCST_3.3G_5Lambda_0.2.csv';
% csvFilePath = 'rcsCST_3.3G_10Lambda.csv';

levelMin = -30;

%% for mesh
thetaDeg = 0:deltaTheta:90;
theta = thetaDeg/180*pi;

phiDeg = 0:deltaPhi:360-deltaPhi;
phi = phiDeg/180*pi;

[theta, phi] = meshgrid(theta, phi);

xMesh = sin(theta).*cos(phi);
yMesh = sin(theta).*sin(phi);
zMesh = cos(theta);

%%
thetaIncident = thetaIncidentDeg/180*pi;

phiIncident = phiIncidentDeg/180*pi;

thetaIncidentMask = zeros(size(thetaDeg));
thetaIncidentMask(find(thetaDeg >= thetaIncidentDeg, 1, 'first')) = 1;

phiIncidentMask = zeros(size(phiDeg));
phiIncidentMask(find(phiDeg >= phiIncidentDeg, 1, 'first')) = 1;

[thetaIncidentMask, phiIncidentMask] = meshgrid(thetaIncidentMask, phiIncidentMask);

IncidentMask = thetaIncidentMask.*phiIncidentMask;

%%
Sa = sinc(a./lambda*(sin(theta).*cos(phi)+sin(thetaIncident)*cos(phiIncident))).*sinc(b./lambda*(sin(theta).*sin(phi)+sin(thetaIncident)*sin(phiIncident)));

C = (1-Gamma)/2;

E_s_theta = C*a*b./lambda*cos(thetaIncident).*cos(theta).*(cos(phiIncident).*sin(phi)-sin(phiIncident).*cos(phi)).*Sa;
E_s_phi = C*a*b./lambda*cos(thetaIncident).*(sin(phiIncident).*sin(phi)+cos(phiIncident).*cos(phi)).*Sa;
E_s_r = 0;

P_s = E_s_theta.^2+E_s_phi.^2+E_s_r^2;

RCS = 10*log10(4*pi*P_s);

rcsFloor = RCS;
rcsFloor(RCS < levelMin) = levelMin;

rcs_x = [RCS(1,end:-1:1) RCS(180/deltaTheta+1,2:end)];
rcs_y = [RCS(90/deltaTheta+1,end:-1:1) RCS(270/deltaTheta+1,2:end)];

%%
csvFile = readmatrix(csvFilePath);

rcsCST = reshape(csvFile(:,3),90/deltaTheta+1,[]);
rcsCST = rcsCST.';
rcsCSTFloor = rcsCST;
rcsCSTFloor(rcsCST < levelMin) = levelMin;

rcsCST_x = [rcsCST(1,end:-1:1) rcsCST(180/deltaTheta+1,2:end)];
rcsCST_y = [rcsCST(90/deltaTheta+1,end:-1:1) rcsCST(270/deltaTheta+1,2:end)];

%%
figure1 = figure;
axes1 = axes('Parent', figure1);

plot(-90:deltaTheta:90, rcs_x, 'LineWidth',1)
hold on
plot(-90:deltaTheta:90, rcsCST_x, 'LineWidth',1,'LineStyle','-.')
hold off
xlim([-90, 90])
ylim([levelMin, -15])
grid on
box on
legend('Proposed','CST')
xlabel('$\theta(^\circ)$','interpreter','latex')
ylabel('RCS(dB/$m^2$)','interpreter','latex')
% title('$\phi=0^\circ$','interpreter','latex')
title('xz-plane','interpreter','latex')


set(axes1, 'GridLineStyle', ':');

% saveas(gcf, 'PatchRCS_CTS.png');

%%
figure2 = figure;
axes2 = axes('Parent', figure2);

plot(-90:deltaTheta:90, rcs_y, 'LineWidth',1)
hold on
plot(-90:deltaTheta:90, rcsCST_y, 'LineWidth',1,'LineStyle','-.')
hold off
ylim([0, 40])
xlim([-90, 90])
ylim([levelMin, -15])
grid on
box on
legend('Proposed','CST')
xlabel('$\theta(^\circ)$','interpreter','latex')
ylabel('RCS(dB/$m^2$)','interpreter','latex')
% title('$\phi=90^\circ$','interpreter','latex')
title('yz-plane','interpreter','latex')

set(axes2, 'GridLineStyle', ':');

% saveas(gcf, 'PatchRCS_CTS_90Deg.png');

%%
figure
surf(xMesh, yMesh, rcsFloor)
shading interp
axis vis3d
axis tight
box on
xlabel('$\sin(\theta) \cos(\phi)$','interpreter','latex')
ylabel('$\sin(\theta) \sin(\phi)$','interpreter','latex')


figure
surf(xMesh, yMesh, rcsCSTFloor)
shading interp
axis vis3d
axis tight
box on
xlabel('$\sin(\theta) \cos(\phi)$','interpreter','latex')
ylabel('$\sin(\theta) \sin(\phi)$','interpreter','latex')


%%
figure
surf((rcsFloor-levelMin).*xMesh, (rcsFloor-levelMin).*yMesh, (rcsFloor-levelMin).*zMesh, (rcsFloor-levelMin))
shading interp
axis vis3d equal
box on
xlim([-15, 15])
ylim([-15, 15])
axis off
view([0, 25])

% saveas(gcf, 'RCS_Proposed.png');
% 
% imageTmp = imread('RCS_Proposed.png');
% imageTmp = imcrop(imageTmp,[910/2-200 656/2-250 400 480]);
% imwrite(imageTmp,'RCS_Proposed.png')

figure
surf((rcsCSTFloor-levelMin).*xMesh, (rcsCSTFloor-levelMin).*yMesh, (rcsCSTFloor-levelMin).*zMesh, (rcsCSTFloor-levelMin))
shading interp
axis vis3d equal
box on
xlim([-15, 15])
ylim([-15, 15])
axis off
view([0, 25])

% saveas(gcf, 'RCS_CST.png');
% 
% imageTmp = imread('RCS_CST.png');
% imageTmp = imcrop(imageTmp,[910/2-200 656/2-250 400 480]);
% imwrite(imageTmp,'RCS_CST.png')