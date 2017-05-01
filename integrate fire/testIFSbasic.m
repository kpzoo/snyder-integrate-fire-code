% Test the IFS code
clearvars 
clc
close all

% Setup data length
runs = 8000;
runsdiv = 1000;

% Folders with data
folder.IFS = 'processed ifs'; % save processed ifs data
folder.SSAload = 'ssa data/gamma2'; % load ssa data
folder.QBload = 'qb data/qb div/first'; % load divided qb data
%folder.QBload = 'qb data/qb div/';
folder.QBsave = 'qb data/qb save'; % save processed qb data
%folder.QBundiv = 'qb data/qb undiv/first'; % qb data to divide

% File name roots for saving generated data or loading existing data
fileRoot.IFS = 'ifs';
fileRoot.loadQB = 'QBdiv';
fileRoot.saveQB = 'QBpcs';

% Reorganise the QB data
%nruns = formatQBdata(runs, runsdiv, folder);

% Run IFS and time its operation
tic;
ifs = runIntegrateFire(folder, fileRoot, runs, runsdiv);
trun = toc/60;
disp(['Single IFS: ' num2str(trun) ' mins']);

% Parameters for simulations
gamma = ifs.gamma;
beta = ifs.beta;

% MSE values
mseEst2 = ifs.mseEst;
mseEst = ifs.mseEst;

% Integrate-fire data
qbArea = ifs.qbArea;
numPest = ifs.numPest;
numPreal = ifs.numPreal;

% Order the beta and corresponding variables
[beta, sortID] = sort(beta);
mseEst = mseEst(sortID);
mseEst2 = mseEst2(sortID);
qbArea = qbArea(sortID);
numPest = numPest(sortID);
numPreal = numPreal(sortID);

% Check gamma
if all(gamma == gamma(1))
    gam = gamma(1);
end

% Plot the MSE across beta at a given gamma
figure;
plot(beta, mseEst)%, beta, mseEst2);
xlabel('\beta');
ylabel('MSE');
title(['Estimates at \gamma = ' num2str(gam)]);

% Plot the photon estimates
figure;
plot(beta, numPreal, 'ko', beta, numPest, 'rs');
xlabel('\beta');
ylabel('photon counts');
title(['Photon counts at \gamma = ' num2str(gam)]);