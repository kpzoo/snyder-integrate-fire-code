% Main script that runs the integrate-fire-Snyder and compared it to the 
% deterministically delayed Synder filter
clearvars
clc
close all

% Assumptions and modifications
% - assumes IPP model
% - compares various filters for the cascade (includes extrinsic noise)
% - assumes a subfolder has Gillespie runs from the light model
% - assumes this data has already been filtered (extrinsic noise only)

%% Simulation parameter setup

% Source directory
thisDir = cd;

% Subfolders of data
gillFol = 'ssa data/gamma10'; % Gillespie source files
extrinFol = 'ipp data'; % Snyder filter at front of cascade

% Get Gillespise SSA source files (only simulations, not filtered)
cd(gillFol);
gillFiles = dir('*.mat');
cd(thisDir);
% No. files to process
flen = length(gillFiles); 

% Options for deterministic delay
tau = 43.3;
limPhoton = 8000; % photon limit


%% Snyder filtering and Integrate-fire

% Variables to store sim parameters
beta = zeros(1, flen);
gamma = zeros(1, flen);
alpha = zeros(1, flen);

% Store MSE
meanEst = zeros(1, flen);
mseEst = zeros(1, flen);
mseEst2 = zeros(1, flen);

% Load each SSA file and process/filter
for i = 1:flen
    % SSA data to be loaded
    gillDataName = gillFiles(i).name;
    cd(gillFol);
    gillData = load(gillDataName, 'outGil', 'params', 'Q');
    cd(thisDir);
    
    % Extract SSA molecular counts and event times and state space
    X = gillData.outGil.X;
    T = gillData.outGil.T;
    Slim = gillData.params.SlimSet;
    Q = gillData.Q;
    
    % Get 3 key simulation settings
    beta(i) = gillData.params.beta;
    gamma(i) = gillData.params.gamma;
    alpha(i) = gillData.params.alpha;
    
    % Remove rest of loaded data
    clearvars('gillData');
    
    % Prepare the SSA data by removing offsets on x2 and T due to transient 
    X(:, 2) = X(:, 2) - X(1, 2);
    T = T - T(1);

    % Limit the photon stream based on limPhoton input
    if X(end, 2) < limPhoton
        warning('Mat:fewPh', ['There are less than ' num2str(limPhoton) ' photons present ']);
    else
        disp(['Clipping Gillespie data to ' num2str(limPhoton) ' photons']);
        idlim = find(X(:, 2) == limPhoton, 1, 'last');
        X = X(1:idlim, :);
        T = T(1:idlim);
    end
    
    % Get event times of observed process (photon times)
    dx2 = [0; diff(X(:, 2))]; % 0 to get correct index
    % Indices of points where new photon comes
    Tevent = T(dx2 == 1);
    
    % Get state space and uniform prior 
    S = diag(Slim.min(1):Slim.max(1));
    lenS = length(S);
    
    % Uniform (discrete) prior and Lam for Snyder solution
    q0 = (1/lenS)*ones(1, lenS);
    Lam = alpha(i)*S;

    % Add deterministic delay and use optimal delayed Snyder filter
    outDet = detDelayProcessing(Tevent, tau, T, X, q0, Q, Lam, lenS, S);
    
    % Get MSE
    meanEst(i) = outDet.meanEst;
    mseEst(i) = outDet.mseEst;
    mseEst2(i) = outDet.mseEst2;
end

%% Process and plot results

% Sort according to beta
[beta, sortID] = sort(beta);
mseEst = mseEst(sortID);
meanEst = meanEst(sortID);
mseEst2 = mseEst2(sortID);

% Get reference of extrinsic noise MSE
cd('ipp data');
ippFile = dir('*.mat');
numIPP = length(ippFile);

% Each file should have same no. beta, also remove the one with 'gamma = 1'
load(ippFile(1).name, 'lenb');

% Declare composite variables
betaIPP = zeros(numIPP, lenb);
mseIPP = zeros(numIPP, lenb);

% Populate variables from every file
for i = 1:numIPP
    IPPstruc = load(ippFile(i).name, 'beta', 'mseComp');
    betaIPP(i, :) =  IPPstruc.beta;
    mseIPP(i, :) =  IPPstruc.mseComp(1, :);
end
cd(thisDir);

% Take average IPP MSE as reference
mseRef = mean(mseIPP);
betaRef = betaIPP(1, :); % all beta should be same

% Check gamma
if all(gamma == gamma(1))
    gam = gamma(1);
end

% Both MSE estimates
figure;
plot(beta, mseEst, '-', beta, mseEst2, 'o');
xlabel('\beta');
ylabel('MSE');
title(['Estimates at \gamma = ' num2str(gam)]);


% Compare to extrinsic noise
figure;
plot(beta, mseEst, betaRef, mseRef);
xlabel('\beta');
ylabel('MSE');
title(['Estimates at \gamma = ' num2str(gam)]);