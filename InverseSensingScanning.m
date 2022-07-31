clear all;close all;clc;

%%
numTrials = 1000;

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
numSnapshotList = 40:10:80;

successList = zeros(length(numSnapshotList), numSnapshotList(end));

for iSnapshot = 1:length(numSnapshotList)

    numSnapshot = numSnapshotList(iSnapshot)

    %%
    numThetaComp = 120;

    %%
    SineThetaComp = -1:2/numThetaComp:1-2/numThetaComp;

    thetaComp = asin(SineThetaComp);
    thetaDegComp = thetaComp/pi*180;

    %%
    anglesObservation = [-30].'/180*pi;

    numObservaton = length(anglesObservation);

    %%
    for numTargets = 1:1:numSnapshot

        numSuccess = 0;

        for iTrial = 1:1:numTrials
            %%
            indexIncident = randperm(numThetaComp, numTargets);
            amplitudeIncident = 10*rand(numTargets, 1);

            %%
            E_IncidentComp = zeros(length(thetaComp),1);
            E_IncidentComp(indexIncident) = amplitudeIncident;


            %%
            SteeringMatrixIncComp = exp(1i*2*pi.*unitsPos.'*sin(thetaComp));

            SteeringVectorDepComp = exp(1i*2*pi.*sin(anglesObservation)*unitsPos);

            %%
            E_Observed = zeros(numSnapshot, 1);

            A = zeros(numSnapshot, numThetaComp);

            for iTest = 1:1:numSnapshot
                currentW = 2*binornd(1, 0.5*ones(1, numUnits))-1;

                %         currentW = exp(1i*2*pi*rand(1,numUnits));

                current_A = C/lambda*exp(-1i*2*pi*r_s)/r_s*SteeringVectorDepComp*diag(currentW)*(SteeringMatrixIncComp);

                A((iTest-1)*numObservaton+1:iTest*numObservaton,:) = current_A;
                E_Observed((iTest-1)*numObservaton+1:iTest*numObservaton) = current_A*E_IncidentComp;
            end

            %%
            Ei_pred = OMP(E_Observed, A, numSnapshot*numObservaton);

            if norm(Ei_pred-E_IncidentComp) < 1e-5
                numSuccess = numSuccess + 1;
            end

        end

        successList(iSnapshot, numTargets) = numSuccess/numTrials;

    end


end
%%
figure
plot(successList(1,:),'Marker','+')
hold on
plot(successList(2,:),'Marker','*')
plot(successList(3,:),'Marker','diamond')
plot(successList(4,:),'Marker','x')
plot(successList(5,:),'Marker','square')
hold off

legend('$T=40$','$T=50$','$T=60$','$T=70$','$T=80$','Interpreter','latex');

xlabel('Sparsity')
ylabel('Probability of exact recovery')
grid on

exportgraphics(gcf, 'InverseSensing.pdf');