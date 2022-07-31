clear all;close all;clc;

%%
C = -2;

lambda = 1;

r_s = 100;

%%
a = 0.1;
b = 0.1;

%%
d = 0.5;

numUnits = 100;
unitsPos = (0:1:numUnits-1)*d;

%%
anglesIncident = [30, 70]/180*pi;
% anglesIncident = [30]/180*pi;
amplitudeIncident = [1.5, 1, 0.5];

anglesSteer = -50/180*pi;

%%
deltaThetaComp = 0.1;

thetaDegComp = -90:deltaThetaComp:90;
thetaComp = thetaDegComp/180*pi;

%%
SteeringVectorComp = exp(1i*2*pi.*sin(thetaComp)*d);

SteeringMatrixIncComp = fliplr(vander(SteeringVectorComp));
SteeringMatrixIncComp = SteeringMatrixIncComp(:,1:numUnits);

SteeringMatrixDepComp = fliplr(vander(SteeringVectorComp));
SteeringMatrixDepComp = SteeringMatrixDepComp(:,1:numUnits);

%%
E_IncidentComp = zeros(length(thetaComp),1);

for iAngle = 1:1:length(anglesIncident)
    currentAngle = anglesIncident(iAngle);
    E_IncidentComp(find(thetaComp >= currentAngle, 1, 'first')) = amplitudeIncident(iAngle);    
end

%%
w = a*b*exp(-1i*2*pi.*sin(anglesIncident(1)).*unitsPos).*exp(-1i*2*pi.*sin(anglesSteer)*unitsPos);

%%
E_ScatteredComp = C/lambda*exp(-1i*2*pi*r_s)/r_s*SteeringMatrixDepComp*diag(w)*(SteeringMatrixIncComp.')*E_IncidentComp;


%% for plotting
deltaThetaPlot = 0.01;

thetaDegPlot = -90:deltaThetaPlot:90;
thetaPlot = thetaDegPlot/180*pi;

%%
SteeringVectorPlot = exp(1i*2*pi.*sin(thetaPlot)*d);

SteeringMatrixIncPlot = fliplr(vander(SteeringVectorPlot));
SteeringMatrixIncPlot = SteeringMatrixIncPlot(:,1:numUnits);

SteeringMatrixDepPlot = fliplr(vander(SteeringVectorPlot));
SteeringMatrixDepPlot = SteeringMatrixDepPlot(:,1:numUnits);

%%
E_IncidentPlot = zeros(length(thetaPlot),1);

for iAngle = 1:1:length(anglesIncident)
    currentAngle = anglesIncident(iAngle);
    E_IncidentPlot(find(thetaPlot >= currentAngle, 1, 'first')) = amplitudeIncident(iAngle);    
end

E_ScatteredPlot = C/lambda*exp(-1i*2*pi*r_s)/r_s*SteeringMatrixDepPlot*diag(w)*(SteeringMatrixIncPlot.')*E_IncidentPlot;


%%
axes1 = axes('Parent',figure);

colororder({'r','b'})

yyaxis left
plot(thetaDegPlot, abs(E_ScatteredPlot),'LineStyle','-','Color',[0.85 0.33 0.1])
ylabel('$| E^s(r^s, \theta) |$','interpreter','latex')

set(gca,'YColor',[0.85 0.33 0.1])
set(axes1,'GridColor',[0.15 0.15 0.15]);

yyaxis right
plot(thetaDegPlot, E_IncidentPlot)
ylabel('$E^i(\theta)$','interpreter','latex')

xlabel('$\theta$','interpreter','latex')
% ylabel('$E^s$','interpreter','latex')
leg1 = legend('Scattered', 'Incident');
set(leg1,'Interpreter','latex');
xlim([-90, 90])
ylim([0, 1.6])

box on
grid on

exportgraphics(gcf, 'ArrayAnomalousReflection.pdf');