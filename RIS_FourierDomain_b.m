clear all;close all;clc;

%%
d = 0.5;

numUnits = 100;
unitsPos = (0:1:numUnits-1)*d;

numThetaPlot = numUnits*50;

%%
anglesIncident = [30, 70]/180*pi;
amplitudeIncident = [1.5, 1];

anglesSteer = [-50]/180*pi;
amplitudeSteer = [3];

%%
SineThetaComp = -1:2/numUnits:1-2/numUnits;

thetaComp = asin(SineThetaComp);
thetaDegComp = thetaComp/pi*180;

%%
SteeringVectorComp = exp(1i*2*pi.*sin(thetaComp)*d);

SteeringMatrixIncComp = fliplr(vander(SteeringVectorComp)).';
SteeringMatrixIncComp = SteeringMatrixIncComp(1:numUnits,:);

SteeringMatrixDepComp = fliplr(vander(SteeringVectorComp));
SteeringMatrixDepComp = SteeringMatrixDepComp(:,1:numUnits);

%%
E_IncidentComp = zeros(length(thetaComp),1);

for iAngle = 1:1:length(anglesIncident)
    currentAngle = anglesIncident(iAngle);
    E_IncidentComp(find(thetaComp >= currentAngle, 1, 'first')) = amplitudeIncident(iAngle);    
end

E_SteerDesired = zeros(length(thetaComp),1);
for iAngle = 1:1:length(anglesSteer)
    currentAngle = anglesSteer(iAngle);
    E_SteerDesired(find(thetaComp >= currentAngle, 1, 'first')) = amplitudeSteer(iAngle);    
end

%%
w = exp(-1i*2*pi.*sin(anglesIncident(1)).*unitsPos).*exp(-1i*2*pi.*sin(anglesSteer)*unitsPos);

w_Inverse = (SteeringMatrixDepComp\E_SteerDesired)./(SteeringMatrixIncComp*E_IncidentComp);

%%
E_ScatteredComp = SteeringMatrixDepComp*diag(w_Inverse)*SteeringMatrixIncComp*E_IncidentComp;


%% for plotting
% deltaThetaPlot = 0.01;
% 
% thetaDegPlot = -90:deltaThetaPlot:90;
% thetaPlot = thetaDegPlot/180*pi;


SineThetaPlot = -1:2/numThetaPlot:1-2/numThetaPlot;

thetaPlot = asin(SineThetaPlot);
thetaDegPlot = thetaPlot/pi*180;


%%
SteeringVectorPlot = exp(1i*2*pi.*sin(thetaPlot)*d);

SteeringMatrixIncPlot = fliplr(vander(SteeringVectorPlot)).';
SteeringMatrixIncPlot = SteeringMatrixIncPlot(1:numUnits,:);

SteeringMatrixDepPlot = fliplr(vander(SteeringVectorPlot));
SteeringMatrixDepPlot = SteeringMatrixDepPlot(:,1:numUnits);

%%
E_IncidentPlot = zeros(length(thetaPlot),1);

for iAngle = 1:1:length(anglesIncident)
    currentAngle = anglesIncident(iAngle);
    E_IncidentPlot(find(thetaPlot >= currentAngle, 1, 'first')) = amplitudeIncident(iAngle);    
end

E_ScatteredPlot = SteeringMatrixDepPlot*diag(w_Inverse)*SteeringMatrixIncPlot*E_IncidentPlot;


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
figure1 = figure;
axes1 = axes('Parent', figure1);
hold(axes1,'on');

plot(thetaDegPlot, E_IncidentPlot)
plot(thetaDegPlot, abs(E_ScatteredPlot))
xlabel('$\theta$','interpreter','latex')
% ylabel('$E^s$','interpreter','latex')
leg1 = legend('Incident', 'Scattered');
set(leg1,'Interpreter','latex');
xlim([-90, 90])
ylim([0, 1.6])

box(axes1,'on');
grid(axes1,'on');
% set(axes1, 'GridLineStyle', ':');

hold(axes1,'off');

saveas(gcf, 'WeightFourierDomain.png');


% set(gca,'units','centimeters')
% pos = get(gca,'Position');
% ti = get(gca,'TightInset');
% 
% set(gcf, 'PaperUnits','centimeters');
% set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
% 
% saveas(gcf, 'WeightFourierDomain.pdf');