clear all;close all;clc;

%%
deltaTheta = 0.02;

%%
d = 0.5;

numUnits = 100;
unitsPos = (0:1:numUnits-1)*d;

%%
anglesIncident = [30, 70]/180*pi;
amplitudeIncident = [1, 2];

anglesSteer = -50/180*pi;

%%
thetaDeg = -90:deltaTheta:90;
theta = thetaDeg/180*pi;

%%

E_Incident = zeros(length(theta),1);

for iAngle = 1:1:length(anglesIncident)
    currentAngle = anglesIncident(iAngle);
    E_Incident(find(theta >= currentAngle, 1, 'first')) = amplitudeIncident(iAngle);    
end

%%
w = exp(-1i*2*pi.*sin(anglesIncident(1)).*unitsPos).*exp(-1i*2*pi.*sin(anglesSteer)*unitsPos);

%%
SteeringVector = exp(1i*2*pi.*sin(theta)*d);

SteeringMatrix_1 = fliplr(vander(SteeringVector)).';
SteeringMatrix_1 = SteeringMatrix_1(1:numUnits,:);

SteeringMatrix_2 = fliplr(vander(SteeringVector));
SteeringMatrix_2 = SteeringMatrix_2(:,1:numUnits);


%%
E_Scattered = 1/numUnits*SteeringMatrix_2*diag(w)*SteeringMatrix_1*E_Incident;



%%
plot(thetaDeg, E_Incident)
hold on
plot(thetaDeg, abs(E_Scattered))
hold off
grid on