clear all;close all;clc;

%%
speedLight = 299792458;

f = 3.3e9;
% f = 0.3e9;

lambda = speedLight/f;

a = 5*lambda;
b = 5*lambda;

Gamma = -1;

thetaIncidentDeg = [15, 45];
phiIncidentDeg = [-45, 135];
amplitudeIncident = [1, 0.5];

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
C = (1-Gamma)/2;

E_s_theta = zeros(size(xMesh));
E_s_phi = zeros(size(xMesh));

for iIncident = 1:1:length(thetaIncidentDeg)
    currentThetaIncident = thetaIncidentDeg(iIncident)/180*pi;
    currentPhiIncident = phiIncidentDeg(iIncident)/180*pi;
    currentAmplitudeIncident = amplitudeIncident(iIncident);
    
    Sa = sinc(a./lambda*(sin(theta).*cos(phi)+sin(currentThetaIncident)*cos(currentPhiIncident))).*sinc(b./lambda*(sin(theta).*sin(phi)+sin(currentThetaIncident)*sin(currentPhiIncident)));

    E_s_theta = E_s_theta+C*a*b./lambda*currentAmplitudeIncident*cos(currentThetaIncident).*cos(theta).*(cos(currentPhiIncident).*sin(phi)-sin(currentPhiIncident).*cos(phi)).*Sa;

    E_s_phi = E_s_phi+C*a*b./lambda*currentAmplitudeIncident*cos(currentThetaIncident).*(sin(currentPhiIncident).*sin(phi)+cos(currentPhiIncident).*cos(phi)).*Sa;

    E_s_r = 0;

end


E_s = sqrt(E_s_theta.^2+E_s_phi.^2+E_s_r^2);

E_s_Normalized = E_s./max(max(E_s));


%%
xSpikes = sin(thetaIncidentDeg/180*pi).*cos(phiIncidentDeg/180*pi);
ySpikes = sin(thetaIncidentDeg/180*pi).*sin(phiIncidentDeg/180*pi);
zSpikes = amplitudeIncident;

%%
figure
surf(xMesh, yMesh, E_s)
shading interp
hold on
stem3(xSpikes,ySpikes,zSpikes,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 0 0],'Marker','square','Color',[1 0 0]);
hold off

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

exportgraphics(gcf, 'PatchResponseMultipleIncidents_3d.pdf');


%%
figure
mesh(E_s.*xMesh, E_s.*yMesh, E_s.*zMesh, E_s)
axis vis3d equal
box on
% axis off
xlim([-0.6, 0.6])
ylim([-0.6, 0.6])
colorbar

exportgraphics(gcf, 'PatchResponseMultipleIncidentsPolar_3d.pdf');