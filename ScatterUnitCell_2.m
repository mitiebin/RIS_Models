clear all;close all;clc;

%%
Gamma = 1;

lambda = 1;
omega = 2*pi/lambda;

%%
a = 0.5;
b = 0.5;

theta_i = 0/180*pi;
phi_i = 270/180*pi;

theta_s = (0:0.01:89)/180*pi;
phi_s = 90/180*pi*ones(size(theta_s));


%%
Sinc_X = sinc(omega*a/pi*sin(theta_s).*cos(phi_s));

Sinc_Y = sinc(omega*b/pi*(sin(theta_s).*sin(phi_s)-sin(theta_i)));

E_s_E_i = 1/2*Gamma*a*b/lambda*sqrt(cos(theta_s).^2.*sin(phi_s).^2+cos(phi_s).^2).*Sinc_X.*Sinc_Y;

RCS = 10*log10(4*pi*E_s_E_i.^2);

%%
figure
plot(0:0.01:89,RCS)

grid on