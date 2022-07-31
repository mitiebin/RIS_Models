clear all;close all;clc;

%%
d = 0.5;

numUnits = 100;
unitsPos = (0:1:numUnits-1)*d;

numThetaComp = 200;
numThetaPlot = 1e4;

%%
numIncidents = 100;

anglesIncident = (rand(numIncidents,1)*180-90)/180*pi;
amplitudeIncident = rand(numIncidents,1);

anglesSteer = -50/180*pi;

%%
% deltaThetaComp = 0.1;
% 
% thetaDegComp = -90:deltaThetaComp:90;
% thetaComp = thetaDegComp/180*pi;
% 
% SteeringVectorComp = exp(1i*2*pi.*sin(thetaComp)*d);

%%
SineThetaComp = -1:2/numThetaComp:1-2/numThetaComp;

thetaComp = asin(SineThetaComp);
thetaDegComp = thetaComp/pi*180;

SteeringVectorComp = exp(1i*2*pi.*sin(thetaComp)*d);

%%
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

%%
wDelay = exp(-1i*2*pi.*sin(anglesIncident(1)).*unitsPos).*exp(-1i*2*pi.*sin(anglesSteer)*unitsPos);

wRandom = sign(rand(numUnits,1)-0.5);


%%
E_ScatteredDelayComp = SteeringMatrixDepComp*diag(wDelay)*SteeringMatrixIncComp*E_IncidentComp;

E_ScatteredRandomComp = SteeringMatrixDepComp*diag(wRandom)*SteeringMatrixIncComp*E_IncidentComp;


%% for plotting
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

E_ScatteredDelayPlot = SteeringMatrixDepPlot*diag(wDelay)*SteeringMatrixIncPlot*E_IncidentPlot;

E_ScatteredRandomPlot = SteeringMatrixDepPlot*diag(wRandom)*SteeringMatrixIncPlot*E_IncidentPlot;


%%
max(abs(E_ScatteredDelayComp))/((sum(abs(E_ScatteredRandomComp)))/length(thetaComp))

max(abs(E_ScatteredDelayPlot))/((sum(abs(E_ScatteredRandomPlot)))/length(thetaPlot))


%%
figure1 = figure;
axes1 = axes('Parent', figure1);
hold(axes1,'on');

plot(thetaDegPlot, abs(E_ScatteredRandomPlot))
plot(thetaDegPlot, abs(E_ScatteredDelayPlot))
plot(thetaDegComp, abs(E_ScatteredRandomComp))

xlabel('$\theta$','interpreter','latex')
% ylabel('$E^s$','interpreter','latex')
leg1 = legend('Incident', 'Scattered');
set(leg1,'Interpreter','latex');
xlim([-90, 90])

box(axes1,'on');
grid(axes1,'on');
set(axes1, 'GridLineStyle', ':');

hold(axes1,'off');
title('$c=0.5$','Interpreter','latex')


% set(gca,'units','centimeters')
% pos = get(gca,'Position');
% ti = get(gca,'TightInset');
% 
% set(gcf, 'PaperUnits','centimeters');
% set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
% 
% saveas(gcf, 'WeightTimeDelay.pdf');