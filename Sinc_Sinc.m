clear all; clc; close all;

%%
theta_i = 70/180*pi;

theta_s = (0:1:90)/180*pi;
phi_s = (0:1:359)/180*pi;


%%
E_i = 1;

r = 1;

lambda = 1;

a = 0.05*lambda;
b = 0.05*lambda;

omega = 2*pi/lambda;

%%
[theta_s, phi_s] = meshgrid(theta_s, phi_s);

%%
Sinc_A = sinc(omega.*a.*sin(theta_s).*cos(phi_s)/(2*pi));
Sinc_B = sinc(omega.*b.*(sin(theta_s).*sin(phi_s)-sin(theta_i))/(2*pi));

Sqrt_A = sqrt(cos(theta_s).^2.*sin(phi_s).^2+cos(phi_s).^2);

E_o = E_i/r*(a*b/lambda).*Sqrt_A.*Sinc_A.*Sinc_B;

%%
figure
mesh(theta_s,phi_s,Sqrt_A)
zlim([0,1])

figure
mesh(theta_s,phi_s,Sinc_A)
zlim([0,1])

figure
mesh(theta_s,phi_s,Sinc_B)
zlim([0,1])