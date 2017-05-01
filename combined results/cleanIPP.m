% File to clean up IPP-Snyder data
clearvars
clc 
close all

% Get the folders of interest
thisDir = cd;
cd('Snyder photon');
fols = dir('batch_*');

% Get basic information (assumed constant across files)
cd(fols(1).name);
% Data files and number - across gamma assumed
files = dir('*.mat');
nFiles = length(files);
% Main info that is assumed constant
info = load(files(1).name, 'gamma', 'lenb', 'k', 'leng');
cd ..

% No. folders 
nFols = length(fols);
leng = info.leng;
gamma = info.gamma;
lenb = info.lenb;

% Combine data for each folder into a single file so 1 per run
for i = 1:nFols
    % Declare variables, ig picks entry in leng for a given file
    betaSet = zeros(leng, lenb);
    mseEst = zeros(leng, lenb);
    meanEst = zeros(leng, lenb);
    gammaSet = zeros(leng, 1);
    alphaSet = zeros(leng, lenb);
    
    % Load folder of data
    cd(fols(i).name);
    % Load data from each file and combine
    for j = 1:nFiles
        tempData = load(files(j).name, 'mseComp', 'meanComp', 'alpha', 'ig', 'beta');
        meanEst(j, :) = tempData.meanComp(1, :);
        mseEst(j, :) = tempData.mseComp(1, :);
        alphaSet(j, :) = tempData.alpha;
        betaSet(j, :) = tempData.beta;
        gammaSet(j) = gamma(tempData.ig);
    end
    cd ..
    % Save file with this data
    save(['completeIPP_' num2str(i)], 'meanEst', 'mseEst', 'alphaSet', 'gammaSet', 'betaSet');
end
cd(thisDir);
        