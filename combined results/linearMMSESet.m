% Calculation of linear MMSE for on-off modulated light
clear all
clc
close all

% Assumptions and modifications
% - calculates steady state linear MMSE
% - works across gamma via k

% Gamma of interest
gamma = [5 10 20 30];
leng = length(gamma);
kset = 1./(100*gamma);

% Set main parameters
beta = 1:0.1:200;
lenb = length(beta);

% Main variables to store
mmseLin = zeros(leng, lenb);

for i = 1:leng
    % Alpha
    k = kset(i);
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
    
    % Store MMSE
    mmseLin(i, :) = mmse;
end

% Remove 'imaginary parts'
mmseLin = abs(mmseLin);
save('linearData.mat', 'mmseLin', 'beta', 'gamma');

% Plot all MMSE curves with gamma
figure;
plot(beta, mmseLin);
xlabel('beta');
ylabel('mmse');
title('Linear MMSE for IPP');
