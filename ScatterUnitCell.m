clear all;close all;clc;

%%
Gamma = 1;

lambda = 1;
omega = 2*pi/lambda;

%%
a = 0.2;
b = 0.2;

theta_i = 0/180*pi;
phi_i = 270/180*pi;

theta_s = (0:1:89)/180*pi;
phi_s = (0:1:359)/180*pi;

[theta_s,phi_s] = meshgrid(theta_s,phi_s);

%%
Sinc_X = sinc(omega*a/pi*sin(theta_s).*cos(phi_s));

Sinc_Y = sinc(omega*b/pi*(sin(theta_s).*sin(phi_s)-sin(theta_i)));


E_s_E_i = 1/2*Gamma*a*b/lambda*sqrt(cos(theta_s).^2.*sin(phi_s).^2+cos(phi_s).^2).*Sinc_X.*Sinc_Y;

RCS = 20*log(E_s_E_i);

%%
figure
mesh(theta_s,phi_s,E_s_E_i)
xlabel('theta')
ylabel('phi')
zlabel('Ratio')


figure
mesh(theta_s,phi_s,RCS)
xlabel('theta')
ylabel('phi')
zlabel('RCS')
