% Code to convert raw QB data into divided form
clearvars
clc
close all

% Assumptions and modifications
% - for several independent replicates with same runs value
% - works with folder structure of createQBsetLoop
% - the same folder sub-structure defined by folQB exists for undiv and div

% Specify the number of QBs and QBs per file desired
runs = 8000;
runsdiv = 1000;

% Define paths to folder
thisDir = cd;
folder.QBundiv = 'qb data/'; % qb data to divide
folder.QBload = 'qb div/'; % save divided qb data here

% Folders names of QBs want to process
cd(folder.QBundiv);
folQB = dir('qb*');
cd(thisDir);
nFiles = length(folQB);

% Confirmation of no. runs
nruns = zeros(1, nFiles);

% Loop across and process QBs in turn
for i = 1:nFiles
    % Specific folders considered
    folTemp.QBundiv = [folder.QBundiv folQB(i).name]; 
    folTemp.QBload = [folder.QBload folQB(i).name];

    % Reorganise the QB data
    nruns(i) = formatQBdata2(runs, runsdiv, folTemp);
    disp(['Formatted ' num2str(i) ' of ' num2str(nFiles)]);
end