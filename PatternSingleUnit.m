clear all;close all;clc;

%%
thetaIncident = 30/180*pi;

bList = [0.1, 0.5, 1, 5];

%%
deltaTheta = 0.01;

thetaDeg = -90:deltaTheta:90;
theta = thetaDeg/180*pi;

%%
Response = zeros(length(bList), length(theta));

for iB = 1:1:length(bList)
    b = bList(iB);
    Response(iB,:) = sinc(b.*(sin(theta)+sin(thetaIncident)));
end


%%
figure1 = figure;
axes1 = axes('Parent', figure1);
plot(thetaDeg, abs(Response(1,:)),'Color',[0 0 1]);
annotation(figure1,'textarrow',[0.75 0.7],[0.875 0.835],'String','$b = \frac{\lambda}{10}$', 'Interpreter','latex','HeadWidth',5,'HeadStyle','vback1','HeadLength',5,'Color',[0 0 1]);
hold(axes1,'on');

plot(thetaDeg, abs(Response(2,:)),'Color',[1 0 0]);
annotation(figure1,'textarrow',[0.65 0.6],[0.71 0.66],'String','$b = \frac{\lambda}{2}$', 'Interpreter','latex','HeadWidth',5,'HeadStyle','vback1','HeadLength',5,'Color',[1 0 0]);

plot(thetaDeg, abs(Response(3,:)),'Color',[0 1 0]);
annotation(figure1,'textarrow',[0.61 0.56],[0.47 0.42],'String','$b = \lambda$', 'Interpreter','latex','HeadWidth',5,'HeadStyle','vback1','HeadLength',5,'Color',[1 0 0]);

plot(thetaDeg, abs(Response(4,:)),'Color',[0 1 0]);
annotation(figure1,'textarrow',[0.29 0.34],[0.40 0.35],'String','$b = 5 \lambda$', 'Interpreter','latex','HeadWidth',5,'HeadStyle','vback1','HeadLength',5);


xlim([-90, 90])
ylim([0, 1.1])
xlabel('$\theta$','interpreter','latex')
% ylabel('$E^s$','interpreter','latex')
% set(leg1,'Interpreter','latex');

box(axes1,'on');
grid(axes1,'on');
set(axes1, 'GridLineStyle', ':');

hold(axes1,'off');
%title('$c=0.5$','Interpreter','latex')


set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

saveas(gcf, 'SincPlot.pdf');