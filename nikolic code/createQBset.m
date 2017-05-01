% Script to batch run the Nikolic model to obtain high numbers of QBs
% Creates the required QBs which must then be time offset based on
% Gillespie simulation of photon times to get the overall summed stream
clearvars
clc
close all

% Assumptions and modifications
% - code will aim to get a max (as stochastic) of totruns QBs

% Set the loop variables
totruns = 10000;
div = 1000; % number of QBs wanted per file (division)
len = round(totruns/div); % number of files
runs = div*ones(1, len); 
% Index for saved files names (data is QBs)
biosimroot = 'QB';

% Run the Nikolic model in a loop and store data
for i = 1:len
    biosimname = [biosimroot num2str(i)];
    % Simulates runs(i) QBs for each i and saves to file name: biosimname
    QBtemp = simQBs(runs(i), biosimname);
    disp('-----------------------------------------------------------------');
    disp(['Completed segment: ' num2str(i)]);
    disp('-----------------------------------------------------------------');
end