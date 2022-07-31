clear all;close all;clc;

%%
b = 0.45;
lambda = 1;


%%
thetaDeg = 160;
theta = thetaDeg/180*pi;


theta_i_deg = -thetaDeg:0.5:0;
theta_ii = theta_i_deg/180*pi;

theta_s = theta + theta_ii;

%%
Angle = cos(theta_ii).*sinc(b/lambda*(sin(theta+theta_ii)+sin(theta_ii)));

f = @(x) pi*b/lambda*cos(x)*sin(theta+2*x)*cos(pi*b/lambda*(sin(theta+x)+sin(x))) - sin(pi*b/lambda*(sin(theta+x)+sin(x))) ;

x = fsolve(f,0)/pi*180

%%
plot(theta_i_deg, abs(Angle))
