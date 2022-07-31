clear all;close all;clc;

%%
d = 0.5;

numUnits = 1000;
unitsPos = (0:1:numUnits-1)*d;

%%
anglesIncident = [30, 70]/180*pi;
amplitudeIncident = [1, 2];

anglesSteer = [-50]/180*pi;
amplitudeSteer = [10];

%%
SineTheta = -1:2/numUnits:1-2/numUnits;

theta = asin(SineTheta);
thetaDeg = theta/pi*180;

%%
SteeringVector = exp(1i*2*pi.*sin(theta)*d);

SteeringMatrix_1 = fliplr(vander(SteeringVector)).';
SteeringMatrix_1 = SteeringMatrix_1(1:numUnits,:);

SteeringMatrix_2 = fliplr(vander(SteeringVector));
SteeringMatrix_2 = SteeringMatrix_2(:,1:numUnits);


%%
E_Incident = zeros(length(theta),1);

for iAngle = 1:1:length(anglesIncident)
    currentAngle = anglesIncident(iAngle);
    E_Incident(find(theta >= currentAngle, 1, 'first')) = amplitudeIncident(iAngle);    
end

E_SteerDesired = zeros(length(theta),1);
for iAngle = 1:1:length(anglesSteer)
    currentAngle = anglesSteer(iAngle);
    E_SteerDesired(find(theta >= currentAngle, 1, 'first')) = amplitudeSteer(iAngle);    
end

%%
w_delay = exp(-1i*2*pi.*sin(anglesIncident(1)).*unitsPos).*exp(-1i*2*pi.*sin(anglesSteer)*unitsPos);

w_Inverse = (SteeringMatrix_2\E_SteerDesired)./(SteeringMatrix_1*E_Incident);


%%
E_Scattered = 1/numUnits*SteeringMatrix_2*diag(w_delay)*SteeringMatrix_1*E_Incident;

E_Scattered2 = SteeringMatrix_2*diag(w_Inverse)*SteeringMatrix_1*E_Incident;

%%
deltaThetaPlot = 0.02;

thetaDegPlot = -90:deltaThetaPlot:90;
thetaPlot = thetaDegPlot/180*pi;

SteeringVectorPlot = exp(1i*2*pi.*sin(thetaPlot)*d);

SteeringMatrix_1_Plot = fliplr(vander(SteeringVectorPlot)).';
SteeringMatrix_1_Plot = SteeringMatrix_1_Plot(1:numUnits,:);

SteeringMatrix_2_Plot = fliplr(vander(SteeringVectorPlot));
SteeringMatrix_2_Plot = SteeringMatrix_2_Plot(:,1:numUnits);


E_IncidentPlot = zeros(length(thetaPlot),1);

for iAngle = 1:1:length(anglesIncident)
    currentAngle = anglesIncident(iAngle);
    E_IncidentPlot(find(thetaPlot >= currentAngle, 1, 'first')) = amplitudeIncident(iAngle);    
end

E_Scattered_Plot = SteeringMatrix_2_Plot*diag(w_Inverse)*SteeringMatrix_1_Plot*E_IncidentPlot;


%%
figure(1)
plot(abs(SteeringMatrix_1*E_Incident))

figure(2)
plot(thetaDeg, E_Incident)
hold on
plot(thetaDeg, abs(E_Scattered2))
hold off
grid on

figure(3)
plot(thetaDegPlot, abs(E_Scattered_Plot))
grid on