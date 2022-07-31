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
numTargets = 10;

numSnapshot = 40;

%%
numThetaComp = 90;


%%
SineThetaComp = -1:2/numThetaComp:1-2/numThetaComp;

thetaComp = asin(SineThetaComp);
thetaDegComp = thetaComp/pi*180;


%%
anglesIncident = randperm(90, numTargets)/180*pi;
amplitudeIncident = 10*rand(1, numTargets);

anglesObservation = [-30].'/180*pi;

numObservaton = length(anglesObservation);

%%
E_IncidentComp = zeros(length(thetaComp),1);

for iAngle = 1:1:length(anglesIncident)
    currentAngle = anglesIncident(iAngle);
    E_IncidentComp(find(thetaComp >= currentAngle, 1, 'first')) = amplitudeIncident(iAngle);
end


%%
SteeringMatrixIncComp = exp(1i*2*pi.*unitsPos.'*sin(thetaComp));

SteeringVectorDepComp = exp(1i*2*pi.*sin(anglesObservation)*unitsPos);

%%
E_Observed = zeros(numSnapshot, 1);

A = zeros(numSnapshot, numThetaComp);

for iTest = 1:1:numSnapshot
    currentW = binornd(1, 0.5*ones(1, numUnits));

    %         currentW = exp(1i*2*pi*rand(1,numUnits));

    current_A = C/lambda*exp(-1i*2*pi*r_s)/r_s*SteeringVectorDepComp*diag(currentW)*(SteeringMatrixIncComp);

    A((iTest-1)*numObservaton+1:iTest*numObservaton,:) = current_A;
    E_Observed((iTest-1)*numObservaton+1:iTest*numObservaton) = current_A*E_IncidentComp;
end

%%
Ei_pred = OMP(E_Observed, A, numSnapshot*numObservaton);

norm(Ei_pred-E_IncidentComp)