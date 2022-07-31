clear all;close all;clc;

%%
speedLight = 299792458;

f = 3.3e9;
% f = 0.3e9;

lambda = speedLight/f;

a = 5*lambda;
b = 5*lambda;

Gamma = -1;

thetaIncidentDeg = 0;
phiIncidentDeg = 0;

%%
deltaTheta = 0.1;
deltaPhi = 0.5;

%% for mesh
thetaDeg = 0:deltaTheta:90;
theta = thetaDeg/180*pi;

phiDeg = 0:deltaPhi:360;
phi = phiDeg/180*pi;

[theta, phi] = meshgrid(theta, phi);

xMesh = sin(theta).*cos(phi);
yMesh = sin(theta).*sin(phi);
zMesh = cos(theta);

%%
thetaIncident = thetaIncidentDeg/180*pi;
phiIncident = phiIncidentDeg/180*pi;
amplitudeIncident = 1;

%%
thetaIncidentMask = zeros(size(thetaDeg));
thetaIncidentMask(find(thetaDeg >= thetaIncidentDeg, 1, 'first')) = 1;

phiIncidentMask = zeros(size(phiDeg));
phiIncidentMask(find(phiDeg >= phiIncidentDeg, 1, 'first')) = 1;

[thetaIncidentMask, phiIncidentMask] = meshgrid(thetaIncidentMask, phiIncidentMask);

IncidentMask = thetaIncidentMask.*phiIncidentMask;


%%
Sa = sinc(a/lambda*(sin(theta).*cos(phi)+sin(thetaIncident)*cos(phiIncident))).*sinc(b/lambda*(sin(theta).*sin(phi)+sin(thetaIncident)*sin(phiIncident)));

C = (1-Gamma)/2;

E_s_theta = C*a*b./lambda*amplitudeIncident*cos(thetaIncident).*cos(theta).*(cos(phiIncident).*sin(phi)-sin(phiIncident).*cos(phi)).*Sa;

E_s_phi = C*a*b./lambda*amplitudeIncident*cos(thetaIncident).*(sin(phiIncident).*sin(phi)+cos(phiIncident).*cos(phi)).*Sa;

E_s_r = 0;

E_s = sqrt(E_s_theta.^2+E_s_phi.^2+E_s_r^2);

E_s_Normalized = E_s./max(max(E_s));

%%
figure
surf(xMesh, yMesh, E_s_Normalized)
hold on
mesh(xMesh, yMesh, IncidentMask)
hold off
shading interp
% axis vis3d
axis tight
xlabel('$\sin(\theta) \cos(\phi)$','interpreter','latex')
ylabel('$\sin(\theta) \sin(\phi)$','interpreter','latex')
box on
% colorbar
% colormap winter

xh = get(gca,'XLabel'); % Handle of the x label
set(xh, 'Units', 'Normalized')
pos = get(xh, 'Position');
set(xh, 'Position',pos.*[1.1,-0.5,1],'Rotation',15)
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[0.7,-0.6,1],'Rotation',-25)

saveas(gcf, 'PatchResponseSingleIncident_3d.png');


set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

% saveas(gcf, 'PatchResponseSingleIncident_3d.pdf');

%%
figure
mesh(E_s.*xMesh, E_s.*yMesh, E_s.*zMesh, E_s)
axis vis3d equal
box on
% xlim([-0.2, 0.2])
% ylim([-0.2, 0.2])