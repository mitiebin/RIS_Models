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
numThetaPlot = numUnits*50;

%%
anglesIncident = [30, 70]/180*pi;
amplitudeIncident = [1.5, 1];

anglesSteer = [-50]/180*pi;
amplitudeSteer = [3]/r_s;

%%
SineThetaComp = -1:2/numUnits:1-2/numUnits;

thetaComp = asin(SineThetaComp);
thetaDegComp = thetaComp/pi*180;

%%
E_SteerDesired = zeros(length(thetaComp),1);
for iAngle = 1:1:length(anglesSteer)
    currentAngle = anglesSteer(iAngle);
    E_SteerDesired(find(thetaComp >= currentAngle, 1, 'first')) = amplitudeSteer(iAngle);    
end


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
beta = numUnits*C/lambda*exp(-1i*2*pi*r_s)/r_s;

w_Inverse = 1/beta*(SteeringMatrixDepComp'*E_SteerDesired)./((SteeringMatrixIncComp.')*E_IncidentComp);

%%
E_ScatteredComp = C/lambda*exp(-1i*2*pi*r_s)/r_s*SteeringMatrixDepComp*diag(w_Inverse)*(SteeringMatrixIncComp.')*E_IncidentComp;


%% for plotting
SineThetaPlot = -1:2/numThetaPlot:1-2/numThetaPlot;

thetaPlot = asin(SineThetaPlot);
thetaDegPlot = thetaPlot/pi*180;


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

E_ScatteredPlot = C/lambda*exp(-1i*2*pi*r_s)/r_s*SteeringMatrixDepPlot*diag(w_Inverse)*(SteeringMatrixIncPlot.')*E_IncidentPlot;


%%
figure(1)
plot(thetaDegPlot, E_IncidentPlot)
hold on
plot(thetaDegComp,abs(E_ScatteredComp),'--gs',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
plot(thetaDegPlot, abs(E_ScatteredPlot))
hold off
grid on

%%
figure(2)
plot(thetaDegComp, E_IncidentComp)
hold on
plot(thetaDegComp,abs(E_ScatteredComp),'--gs',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
hold off
grid on


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

exportgraphics(gcf, 'ArrayFourierDomain.pdf');