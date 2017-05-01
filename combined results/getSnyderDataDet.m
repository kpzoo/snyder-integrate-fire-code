% Function to extract Snyder data
function [betaSny, gammaSny, meanSny, mseSny, L, U, mseRef, betaRef] = getSnyderDataDet(locFol, thisDir, numIPP, lenb, leng)

% Data location
cd(locFol);
ippFile = dir('*.mat');

% Declare composite variables, sorted for gamma
betaSny = zeros(numIPP, lenb, leng);
gammaSny = zeros(numIPP, lenb, leng);
meanSny = zeros(numIPP, lenb, leng);
mseSny = zeros(numIPP, lenb, leng);

% Populate variables from every file
for i = 1:numIPP
    IPPstruc = load(ippFile(i).name, 'mseEst', 'meanEst', 'beta', 'gamma');
    for j = 1:leng
        betaSny(i, :, j) =  IPPstruc.beta(j, :);
        meanSny(i, :, j) =  IPPstruc.meanEst(j, :);
        mseSny(i, :, j) =  IPPstruc.mseEst(j, :);
        gammaSny(i, :, j) =  IPPstruc.gamma(j, :);
    end
end
cd(thisDir);

% Take average IPP MSE as reference for each gamma
mseRef = zeros(leng, lenb);
L = zeros(leng, lenb);
U = zeros(leng, lenb);
for i = 1:leng
    mseRef(i, :) = mean(mseSny(:, :, i));
    % Lower and upper bounds
    L(i, :) = mseRef(i, :) - min(mseSny(:, :, i));
    U(i, :) = max(mseSny(:, :, i)) - mseRef(i, :);
end
betaRef = betaSny(1, :, 1); % all beta should be same

% Go home
cd(thisDir);