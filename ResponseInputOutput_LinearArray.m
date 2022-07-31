clear all;close all;clc;

%%
phiDelta = 0.2;

numElements = 20;

distance = 0.49;

phiArr = 60/180*pi;
phiDep = 60/180*pi;

%%
locationElements = (0:1:numElements-1)*distance;

Weight = exp(-1i*2*pi.*locationElements.*(sin(phiArr)+sin(phiDep))).';

% Weight = exp(1i*2*pi*50*randn(size(locationElements))).';
% Weight = ones(size(locationElements)).';

%%
phiIncomingVector = (-90:phiDelta:90)./180*pi;
phiOutgoingVector = (-90:phiDelta:90)./180*pi;

[phiIncomingMesh, phiOutgoingMesh] = meshgrid(phiIncomingVector, phiOutgoingVector);

phiIncoming = reshape(phiIncomingMesh, [], 1);
phiOutgoing = reshape(phiOutgoingMesh, [], 1);

distance_sin_phiIncoming = kron(locationElements, sin(phiIncoming));
distance_sin_phiOutgoing = kron(locationElements, sin(phiOutgoing));

%%
Response = exp(1i*2*pi.*(distance_sin_phiIncoming+distance_sin_phiOutgoing))*Weight;
Response = reshape(Response,size(phiIncomingMesh))./numElements;

%%
[phiIncomingPlot, phiOutgoingPlot] = meshgrid(-90:phiDelta:90, -90:phiDelta:90);


figure(1)
mesh(phiIncomingPlot,phiOutgoingPlot,abs(Response))
xlabel('$\phi_i$ (degree)','interpreter','latex')
ylabel('$\phi_s$ (degree)','interpreter','latex')
zlabel('$f(\phi_i,\phi_s)$','interpreter','latex')
axis tight
rotate3d on
% axis equal
view(45, 30);
% view(0, 90);
box on

figure(2)
mesh(sin(phiIncomingVector),sin(phiOutgoingVector),abs(Response))
xlabel('$\sin(\phi_i)$','interpreter','latex')
ylabel('$\sin(\phi_s)$','interpreter','latex')