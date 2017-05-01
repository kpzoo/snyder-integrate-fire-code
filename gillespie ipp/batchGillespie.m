% Script to runs Gillespie simulations across different alpha for a given
% gamma of interest (beta range set)
clearvars
clc
close all

% Assumptions and modifications
% can do multi-state MC simulations, and 2 state IPP

% Define range of beta and a fixed gamma and obtain rate parameters
%beta = [1 10:10:200];
beta = 30;
gamma = 1;
k = 1/(100*gamma);
alpha = beta*k;
lenb = length(beta);
simroot = 'IPP';

% Output parameters from Gillespie SSA
paramSet = cell(1, lenb);

% Maximum state and molecular count jump (1 => nearest neighbour)
stateMax = 1;
jumpMax = 1;

% Loop across settings
for ib = 1:lenb
    % Set name for saved data
    simname = [simroot num2str(ib)];
    
    % Main code for running SSA
    paramSet{ib} = runGillespie(alpha(ib), k, simname, stateMax, jumpMax);
    
    % Display progress
    disp(['Completed run: ' num2str(ib) ' of ' num2str(lenb)]);
end

% Save data 
save(['batch_' num2str(gamma) '.mat']);

% Variables to store MSE (for testing so temporary)
mseComp = zeros(2, lenb);
meanComp = zeros(2, lenb);

% Process results
for ib = 1:lenb
    % Get MSE values
    mseComp(1, ib) = paramSet{ib}.x1stats1.vals(3);
    mseComp(2, ib) = paramSet{ib}.x1stats2.vals(3);
    meanComp(1, ib) = paramSet{ib}.x1stats1.vals(1);
    meanComp(2, ib) = paramSet{ib}.x1stats2.vals(1);
end


