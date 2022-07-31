clear all;close all;clc;

%%
deltaTheta = 0.5;

d = 0.5;

numUnits = 20;

%%
thetaDeg = -90:deltaTheta:90;
theta = thetaDeg/180*pi;

%%
anglesIncident = [30, 50]/180*pi;

E_Incident = zeros(length(theta),1);

for iAngle = 1:1:length(anglesIncident)
    currentAngle = anglesIncident(iAngle);
    E_Incident(find(theta > currentAngle, 1, 'first')) = 1;    
end

%%
anglesSteer = [-20, -30]/180*pi;

E_Steer = zeros(length(theta),1);


for iAngle = 1:1:length(anglesSteer)
    currentAngle = anglesSteer(iAngle);
    E_Steer(find(theta > currentAngle, 1, 'first')) = 1;    
end

indexSteer = find(E_Steer);


%%
SteeringVector = exp(1i*2*pi.*sin(theta)*d);

SteeringMatrix_1 = fliplr(vander(SteeringVector)).';
SteeringMatrix_1 = SteeringMatrix_1(1:numUnits,:);

SteeringMatrix_2 = fliplr(vander(SteeringVector));
SteeringMatrix_2 = SteeringMatrix_2(:,1:numUnits);


%%
cvx_begin
    variable w(numUnits) complex
    minimize( norm( w , 2 ) )
    subject to
        E_Scattered = SteeringMatrix_2*diag(w)*SteeringMatrix_1*E_Incident;
        E_Scattered(indexSteer) == ones(length(indexSteer),1);
cvx_end

%%
plot(thetaDeg, abs(E_Scattered))
hold on
plot(thetaDeg, E_Steer)
grid on