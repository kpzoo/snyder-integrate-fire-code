% Process the IFS results, 1 QB run per SSA
clearvars
clc
close all

% Assumptions and modifications
% - all files assumed to have same beta and gamma
% - linear MMSE generated 
% - Snyder optimised for deterministic delay in a folder
% - IFS on deterministic Nikolic model in a folder
% - IFS on stochastic Nikolic model in a folder
% - Snyder on just photon noise in a folder

% Home
thisDir = cd;

%% Snyder with only photon noise (optimal)

% Load main IPP data
cd('ipp data');
ippFile = dir('*.mat');
numIPP = length(ippFile);

% Each file should have same beta and gamma compilations
load(ippFile(1).name, 'gammaSet', 'betaSet');
leng = length(gammaSet);
lenb = length(betaSet);
cd(thisDir);

% Extract all data and organise for plotting
[betaIPP, gammaIPP, meanIPP, mseIPP, Lipp, Uipp, mseRef, betaRef] = ...
    getSnyderData('ipp data', thisDir, numIPP, lenb, leng);
    
%% IFS full noise data handling - 1 mat file per QB and SSA run

% Extract and restructure all data
[betaIFS, qbArea, meanIFS, mseIFS, Lifs, Uifs, mseIFSRef, qbRef, Lqb, Uqb] = ...
    getIFSData('ifs data', thisDir, numIPP, leng, lenb);

%% IFS noiseless data handling - 1 mat file per QB and SSA run

% Extract and restructure all data
[betaIFSNo, qbAreaNo, meanIFSNo, mseIFSNo, LifsNo, UifsNo, mseIFSRefNo, qbRefNo,...
    LqbNo, UqbNo] = getIFSData('no noise data', thisDir, numIPP, leng, lenb);

%% Deterministic delay handling - 2 mat files

% Extract all data and organise for plotting
[betaDet, gammaDet, meanDet, mseDet, Ldet, Udet, mseDetRef, betaDetRef] = ...
    getSnyderDataDet('det data', thisDir, numIPP, lenb, leng);


%% Linear MMSE

% Load linear data, same for all gamma
linData = load('linearData.mat');
betaLin = linData.beta;
mmseLin = linData.mmseLin;

%% Plotting of all results

% Define gammaRef
gammaRef = zeros(1, leng);
for i = 1:leng
    gammaRef(i) = unique(gammaIPP(:, :, i));
end

% MSE trajectory error of both IFS models with bounds
eIFS = mseIFS - mseIFSNo;
% Take average IPP MSE as reference for each gamma
eRef = zeros(leng, lenb);
Le = zeros(leng, lenb);
Ue = zeros(leng, lenb);
for i = 1:leng
    eRef(i, :) = mean(eIFS(:, :, i));
    % Lower and upper bounds
    Le(i, :) = eRef(i, :) - min(eIFS(:, :, i));
    Ue(i, :) = max(eIFS(:, :, i)) - eRef(i, :);
end

% MSE difference between stochastic and deterministic Nikolic
figure;
for i = 1:leng
    subplot(ceil(leng/2), 2, i);
    % Errorbar with errors given by [L U]
    errorbar(betaRef, eRef(i, :), Le(i, :), Ue(i, :));
    box off
    xlabel(['\beta, \gamma = ' num2str(gammaRef(i))]);
    ylabel('MSE difference');
    xlim([betaRef(1) betaRef(end)])
end

% QB average area plots for both models
figure;
for i = 1:leng
    subplot(ceil(leng/2), 2, i);
    % Errorbar with errors given by [L U]
    errorbar(betaRef, qbRef(i, :), Lqb(i, :), Uqb(i, :), 'linewidth', 2);
    hold on
    errorbar(betaRef, qbRefNo(i, :), LqbNo(i, :), UqbNo(i, :), 'linewidth', 2);
    hold off
    box off
    xlabel(['\beta, \gamma = ' num2str(gammaRef(i))]);
    ylabel('QB mean area');
    xlim([betaRef(1) betaRef(end)])
end

% Do linear MMSE, IFS plots
figure;
plot(betaRef, mean(mseRef), 'k', 'linewidth', 2);
hold on
plot(betaLin, mmseLin(1, :), 'b', 'linewidth', 2); % all mmse curves same
for i = 1:leng
    hold on
    h = errorbar(betaRef, mseIFSRef(i, :), Lifs(i, :), Uifs(i, :), '--', 'linewidth', 2);
    set(h, 'color', [0.8 0.8 0.8]);
    %plot(betaRef, mseIFSRef(i, :), 'r', 'linewidth', 2);
end
hold off
xlim([betaRef(1) betaRef(end)]);
xlabel('\beta');
ylabel('MSE');
legend('non-linear MMSE (photon noise)', 'linear MMSE (photon noise)',...
    'Integrate-Fire-Snyder (photon and cascade noise)', 'location', 'best');
box off

% Main plot - combines all data except linear MMSE
figure;
for i = 1:leng
    subplot(ceil(leng/2), 2, i);
    % Errorbar with errors given by [L U]
    errorbar(betaRef, mseRef(i, :), Lipp(i, :), Uipp(i, :), 'linewidth', 2);
    hold on
    errorbar(betaRef, mseIFSRef(i, :), Lifs(i, :), Uifs(i, :), 'linewidth', 2);
    errorbar(betaRef, mseIFSRefNo(i, :), LifsNo(i, :), UifsNo(i, :), 'linewidth', 2);
    errorbar(betaRef, mseDetRef(i, :), Ldet(i, :), Udet(i, :), 'linewidth', 2);
    hold off
    box off
    xlabel(['\beta (\gamma = ' num2str(gammaRef(i)) ')']);
    ylabel('MSE');
    legend('Snyder, photon noise', 'Integrate-fire, all noise', 'Integrate-fire, delay', 'Snyder, delay');
    xlim([betaRef(1) betaRef(end)])
end


% Save final results with ID of data/time
tdate = datetime('now');
save('combinedRepeat.mat');