clear all;close all;clc;

%%
d = 0.5;

numUnits = 100;
unitsPos = (0:1:numUnits-1)*d;

numThetaComp = 2e3;

numTrials = 1000;

%%
anglesIncident = [30]/180*pi;
amplitudeIncident = [1, 0.5];

anglesSteer = -50/180*pi;

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
E_Scattered_100Trials = zeros(numThetaComp, numTrials);

for iTrial = 1:1:numTrials
     iTrial
%     wRandom = ones(numUnits,1);
%     wRandom(randsample(numUnits,numUnits/2)) = -1;
    wRandom = sign(rand(numUnits,1)-0.5);
%     sum(wRandom)
    E_Scattered_100Trials(:,iTrial) = SteeringMatrixDepComp*diag(wRandom)*SteeringMatrixIncComp*E_IncidentComp;
end

Abs_E_Mean = sqrt(sum(abs(E_Scattered_100Trials).^2, 2)/numTrials);

%%
wDelay = exp(-1i*2*pi.*sin(anglesIncident(1)).*unitsPos).*exp(-1i*2*pi.*sin(anglesSteer)*unitsPos);

E_ScatteredDelayComp = SteeringMatrixDepComp*diag(wDelay)*SteeringMatrixIncComp*E_IncidentComp;



%%
figure1 = figure;
axes1 = axes('Parent', figure1);


plot(thetaDegComp, Abs_E_Mean)
hold(axes1,'on');
plot(thetaDegComp, abs(E_ScatteredDelayComp))

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