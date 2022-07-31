clear all;close all;clc;

%%
speedLight = 299792458;

f = 3.3e9;
% f = 0.3e9;

lambda = speedLight/f;

a = 5*lambda;
b = 5*lambda;

thetaIncidentDeg = 0;
phiIncidentDeg = 0;

Gamma = -1;
C = (1-Gamma)/2;

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

thetaIncidentMask = zeros(size(thetaDeg));
thetaIncidentMask(find(thetaDeg >= thetaIncidentDeg, 1, 'first')) = 1;

phiIncidentMask = zeros(size(phiDeg));
phiIncidentMask(find(phiDeg >= phiIncidentDeg, 1, 'first')) = 1;

[thetaIncidentMask, phiIncidentMask] = meshgrid(thetaIncidentMask, phiIncidentMask);

IncidentMask = thetaIncidentMask.*phiIncidentMask;


%%
Sa = sinc(a./lambda*(sin(theta).*cos(phi)+sin(thetaIncident)*cos(phiIncident))).*sinc(b./lambda*(sin(theta).*sin(phi)+sin(thetaIncident)*sin(phiIncident)));

E_s_theta = C*a*b./lambda*cos(thetaIncident).*cos(theta).*(cos(phiIncident).*sin(phi)-sin(phiIncident).*cos(phi)).*Sa;

E_s_phi = C*a*b./lambda*cos(thetaIncident).*(sin(phiIncident).*sin(phi)+cos(phiIncident).*cos(phi)).*Sa;

E_s_r = 0;

P_s = E_s_theta.^2+E_s_phi.^2+E_s_r^2;

RCS = 10*log10(4*pi*P_s);

%%
figure
surf(xMesh, yMesh, RCS)

shading interp
% axis vis3d
% axis tight
xlabel('$\sin(\theta) \cos(\phi)$','interpreter','latex')
ylabel('$\sin(\theta) \sin(\phi)$','interpreter','latex')

zlim([-50 max(max(RCS))])

% box on
% colorbar
% colormap winter

% saveas(gcf, 'PatchResponseSingle_3d.pdf');
