% Script to runs Gillespie simulations across different alpha for a given
% gamma of interest (beta range set)
clearvars
clc
close all

% Assumptions and modifications
% - now does additional repeats and includes folder structure
% - added skip_Sny
% - modded version of batchGillespie for range of gamma as well
% - can do multi-state MC simulations, and 2 state IPP

% No. independent runs
nIndep = 5;
thisDir = cd;

for ind = 1:nIndep
    
    % Folder to store data
    folStore = ['ssa_' num2str(ind)];
    batStore = ['batch_' num2str(ind)];
    
    % Define range of beta and a fixed gamma and obtain rate parameters
    beta = [1 10:10:200];
    gamma = [5 10 20 30];
    k = 1./(100*gamma);
    lenb = length(beta);
    leng = length(gamma);
    simroot = ['IPP_' num2str(ind) '_'];
    
    % Maximum state and molecular count jump (1 => nearest neighbour)
    stateMax = 1;
    jumpMax = 1;
    
    % Skip Snyder filtering
    skip_Sny = 1;
    disp('No Snyder filtering');
    
    for ig = 1:leng
        % Output parameters from Gillespie SSA
        paramSet = cell(1, lenb);
        
        % Variables to store MSE (for testing so temporary)
        mseComp = zeros(2, lenb);
        meanComp = zeros(2, lenb);
        
        % Values of alpha
        alpha = beta*k(ig);
        
        % Loop across settings
        for ib = 1:lenb
            % Set name for saved data 
            simname = [simroot num2str(ig) '_' num2str(ib)];
            
            % Main code for running SSA
            paramSet{ib} = runGillespie2(alpha(ib), k(ig), simname, stateMax, jumpMax, skip_Sny, folStore);
            
            % Get MSE values
            if ~skip_Sny
                mseComp(1, ib) = paramSet{ib}.x1stats1.vals(3);
                mseComp(2, ib) = paramSet{ib}.x1stats2.vals(3);
                meanComp(1, ib) = paramSet{ib}.x1stats1.vals(1);
                meanComp(2, ib) = paramSet{ib}.x1stats2.vals(1);
            end
            
            % Display progress
            disp(['Completed run: ' num2str(ib) ' of ' num2str(lenb)]);
        end
        
        % Save data
        if exist(batStore, 'dir') ~= 7
            % Folder does not exist unless 7 returned
            mkdir(batStore);
        end
        cd(batStore);
        save(['batch_' num2str(ig) '.mat']);
        cd(thisDir);
    end
end