% Function that loads and extracts IFS data
function [betaIFS, qbArea, meanIFS, mseIFS, Lifs, Uifs, mseIFSRef, qbRef, Lqb, Uqb] =...
    getIFSData(locFol, thisDir, numIPP, leng, lenb)

% Load main IFS data
cd(locFol);
ifsFile = dir('*.mat');
lenIFS = length(ifsFile);

% Quick check of gamma and runs
ifsFol = load(ifsFile(1).name, 'QBfol', 'SSAfol');
if length(ifsFol.QBfol) ~= numIPP || length(ifsFol.SSAfol) ~= leng
    error('Simulation lengths inconsistent');
end

% Declare composite variables, sorted for gamma
betaIFS = zeros(numIPP, lenb, leng);
qbArea = zeros(numIPP, lenb, leng);
meanIFS = zeros(numIPP, lenb, leng);
mseIFS = zeros(numIPP, lenb, leng);
gammaIFS = zeros(numIPP, lenb, leng);

% Populate variables from every file
for i = 1:lenIFS
    IFSstruc = load(ifsFile(i).name, 'mseEst', 'meanEst', 'beta', 'qbArea', 'gamma');
    for j = 1:leng
        betaIFS(i, :, j) =  IFSstruc.beta{j};
        meanIFS(i, :, j) =  IFSstruc.meanEst{j};
        mseIFS(i, :, j) =  IFSstruc.mseEst{j};
        qbArea(i, :, j) =  IFSstruc.qbArea{j};
        gammaIFS(i, :, j) =  IFSstruc.gamma{j};
    end
end
cd(thisDir);

% Take average IPP MSE as reference for each gamma
mseIFSRef = zeros(leng, lenb);
Lifs = zeros(leng, lenb);
Uifs = zeros(leng, lenb);
for i = 1:leng
    mseIFSRef(i, :) = mean(mseIFS(:, :, i));
    % Lower and upper bounds
    Lifs(i, :) = mseIFSRef(i, :) - min(mseIFS(:, :, i));
    Uifs(i, :) = max(mseIFS(:, :, i)) - mseIFSRef(i, :);
end

% Do average and bounds for area
% Take average IPP MSE as reference for each gamma
qbRef = zeros(leng, lenb);
Lqb = zeros(leng, lenb);
Uqb = zeros(leng, lenb);
for i = 1:leng
    qbRef(i, :) = mean(qbArea(:, :, i));
    % Lower and upper bounds
    Lqb(i, :) = qbRef(i, :) - min(qbArea(:, :, i));
    Uqb(i, :) = max(qbArea(:, :, i)) - qbRef(i, :);
end