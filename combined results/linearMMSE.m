% Included calculation of steasy state mmse

% Calculation of linear MMSE for on-off modulated light
clear all
clc
close all

% Set main parameters
k = 1/2000;
beta = 1:0.1:200;
alpha = k*beta;

% Initialise storage variables
y = zeros(size(alpha));
mmse_inf = y;
t = 10000;

% Loop across parameters and calculate trajectory of covariance function y
for j = 1:length(alpha)
    a = 2*alpha(j)*k;
    b = 4*k;
    c = 1 + (b^2)/(4*a);
    d = b/(2*a);
    G = atanh(sqrt(a/c)*(1/(4*k) + d));
    y(j) = sqrt(c/a)*tanh(sqrt(c/a)*a*t + G) - d;
    mmse_inf(j) = k*(sqrt(c/a) - d);
end

% Obtain normalised mmse and the mmse for the rate process
mmse = k*y;
mmse_lam = (alpha.*alpha).*mmse;

% Plot results across beta
figure;
plot(beta, mmse);
xlabel('beta');
ylabel('mmse');
title('MMSE obtained from optimal linear filter for symmetric 2 state MC');