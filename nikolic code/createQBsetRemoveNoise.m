% Script to batch run the Nikolic model to obtain high numbers of QBs
% Creates the required QBs which must then be time offset based on
% Gillespie simulation of photon times to get the overall summed stream
clearvars
clc
close all

% Assumptions and modifications
% - set flags to remove noise (do deterministic calculations)
% - code will aim to get a max (as stochastic) of totruns QBs

tic;

% Set the loop variables
totruns = 10000;
div = 1000; % number of QBs wanted per file (division)
len = round(totruns/div); % number of files
runs = div*ones(1, len); 

% % Flags control noise [flagM flagG flagP flagD flagT]
flags = [0 0 0 0 0]; % all deterministic

% Index for saved files names (data is QBs)
flagStr = num2str(flags);
flagStr = flagStr(~isspace(flagStr));
biosimroot = ['QB' flagStr '_'];

% Run the Nikolic model in a loop and store data
for i = 1:len
    biosimname = [biosimroot num2str(i)];
    % Simulates runs(i) QBs for each i and saves to file name: biosimname
    QBtemp = simQBRemoveNoise(runs(i), biosimname, flags);
    disp('-----------------------------------------------------------------');
    disp(['Completed segment: ' num2str(i)]);
    disp('-----------------------------------------------------------------');
end

% Time the runs
tqbSim = toc/60;
disp(['Sumulated ' num2str(totruns) ' QBs in ' num2str(tqbSim) ' mins']);