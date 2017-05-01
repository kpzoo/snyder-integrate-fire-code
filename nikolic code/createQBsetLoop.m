% Script to batch run the Nikolic model to obtain high numbers of QBs
% Creates the required QBs which must then be time offset based on
% Gillespie simulation of photon times to get the overall summed stream
clearvars
clc
close all

% Assumptions and modifications
% - runs a loop and assigns data to made folders
% - code will aim to get a max (as stochastic) of totruns QBs

% Looping of the runs (replicates)
nRep = 7;
thisDir = cd;

% Noise setting
remNoise = 1;
if remNoise
    disp('Simulating noiseless QBs');
end
flags = [0 0 0 0 0];

% Set the loop variables
totruns = 10000;
div = 1000; % number of QBs wanted per file (division)
len = round(totruns/div); % number of files
runs = div*ones(1, len); 
% Index for saved files names (data is QBs)
biosimroot = 'QB';
locRoot = 'qb data';

% Create every needed folder
for j = 1:nRep
    % Define save location and create folder there if it does not exist
    folStore = ['qb' num2str(j)];
    cd(locRoot);
    if exist(folStore, 'dir') ~= 7
        % Folder does not exist unless 7 returned
        mkdir(folStore);
    end
    cd(thisDir);
end

parfor j = 1:nRep
    % Folder to save data 
    folStore = ['qb' num2str(j)];
    saveLoc = [locRoot '/' folStore];
    
    % Run the Nikolic model in a loop and store data
    for i = 1:len
        biosimname = [biosimroot num2str(i)];
        
        % Simulates runs(i) QBs for each i and saves to file name: biosimname
        if ~remNoise
            QBtemp = simQBs2(div, biosimname, saveLoc);
        else
            QBtemp = simQBRemoveNoise2(div, biosimname, flags, saveLoc);
        end
        
        disp(['Completed segment: ' num2str(i)]);
    end
    disp('-----------------------------------------------------------------');
    disp(['Completed loop: ' num2str(j)]);
    disp('-----------------------------------------------------------------');
end