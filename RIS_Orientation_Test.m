clear all;close all;clc;

%%
b = 0.45;
lambda = 1;

%%
thetaStep = 0.2;
thetaDeg = 1:thetaStep:179;

theta_i_opt = zeros(size(thetaDeg));

for iTheta = 1:1:length(thetaDeg)
    theta = thetaDeg(iTheta)/180*pi;
    f = @(x) pi*b/lambda*cos(x)*sin(theta+2*x)*cos(pi*b/lambda*(sin(theta+x)+sin(x))) - sin(pi*b/lambda*(sin(theta+x)+sin(x)));
    theta_i_opt(iTheta) = fsolve(f,0)/pi*180;
end

%%
theta_i_opt_plot = zeros(size(thetaDeg));

for iTheta = 1:1:length(thetaDeg)
    currentThetaDeg = thetaDeg(iTheta);
    theta = currentThetaDeg/180*pi;
    
    theta_i_deg = -currentThetaDeg:1:0;
    theta_i = theta_i_deg/180*pi;

    theta_s = theta + theta_i;
    
    Coeff = abs(cos(theta_i).*sinc(b/lambda*(sin(theta_s)+sin(theta_i))));
    
    index = find(Coeff==max(Coeff),1,'first');
    
    theta_i_opt_plot(iTheta) = theta_i_deg(index);
end
%%
plot(thetaDeg,theta_i_opt)


%%
diff_theta_i_opt = abs(theta_i_opt-theta_i_opt_plot);
max(diff_theta_i_opt)