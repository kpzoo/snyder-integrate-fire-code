% Code to convert raw QB data into divided form
clearvars
clc
close all

% Assumptions and modifications
% - the same folder sub-structure defined by folQB exists for undiv and div

% Specify the number of QBs and QBs per file desired
runs = 8000;
runsdiv = 1000;

% Folders names of QBs want to process
folQB = {'first', 'second', 'third'};
nFiles = length(folQB);

% Define paths to folder
folder.QBundiv = 'qb data/'; % qb data to divide
folder.QBload = 'qb div/'; % save divided qb data here

% Confirmation of no. runs
nruns = zeros(1, nFiles);

% Loop across and process QBs in turn
for i = 1:nFiles
    % Specific folders considered
    folTemp.QBundiv = [folder.QBundiv folQB{i}]; 
    folTemp.QBload = [folder.QBload folQB{i}];

    % Reorganise the QB data
    nruns(i) = formatQBdata(runs, runsdiv, folTemp);
end