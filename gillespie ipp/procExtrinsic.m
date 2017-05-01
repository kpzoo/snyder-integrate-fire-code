% Script to extract data from batch files produced from Gillespie sims
clearvars
clc
close all

% Assumptions and modifications
% - produces plot of intrinsic noise against beta and gamma
% - folder with files batch_i, i indicates gamma entry
% - beta consistent across files

% Subfolder name and code source directory
fol = 'batch extrinsic';
thisDir = cd;

% Get number of files in folder
cd(fol);
files = dir('*.mat');
cd(thisDir);
nFiles = length(files);

% Desired variables from each file
varList = {'gamma', 'beta', 'ig', 'mseComp', 'meanComp', 'k', 'ib'};

% Variables of interest
gammaIn = zeros(1, nFiles);
betaIn = cell(1, nFiles);
kIn = zeros(1, nFiles);
mseIn = cell(1, nFiles);
mIn = cell(1, nFiles);
numBeta = zeros(1, nFiles);

% Loop across files getting gamma value and MSE
for i = 1:nFiles
    cd(fol);
    % Load variables from simulations
    load(files(i).name, varList{:});
    cd(thisDir);
    % Assign values for sim parameters
    gammaIn(i) = gamma(ig);
    betaIn{i} = beta;
    kIn(i) = k(ig);
    numBeta(i) = ib;
    % Take MSE and mean from first MSE calculation
    mseIn{i} = mseComp(1, :);
    mIn{i} = meanComp(1, :);
    % Clean up workspace
    clear(varList{:});
    disp(['Finished ' num2str(i) ' of ' num2str(nFiles)]);
end

% Reshape cell data into matrices, assuming beta length same
if ~all(numBeta == numBeta(1))
    error('Sims not batched across the same beta');
else
    % Data across fixed beta set
    nB = numBeta(1);
    % Convert cells and reshape
    beta = cell2mat(betaIn);
    beta = reshape(beta, [nB, nFiles]);
    mse = cell2mat(mseIn);
    mse = reshape(mse, [nB, nFiles]);
    m = cell2mat(mIn);
    m = reshape(m, [nB, nFiles]);
end

% Plot MSE and mean across gamma and beta
figure;
plot(beta, mse, 'linewidth', 2);
box off
xlabel('\beta');
ylabel('mse');
title('Extrinsic noise in model (5 runs)');

% Save figure and data
saveas(gcf, 'extrinsic');
save('procExtr.mat');