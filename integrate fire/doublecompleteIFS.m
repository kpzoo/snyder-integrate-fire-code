% Run IFS code on several repeat photon streams at each gamma
clearvars
clc
close all

% Assumptions and modifications
% - runs across independent ssa data now
% - assume the QBs have already been divided (qb div is populated)
% - assume a folder structure for repeat photon sets
% - loop 1 is on gamma, 2 is on photons and 3 on beta

%% Setup the problem and folder structure

tic;

% Setup data length
runs = 8000;
runsdiv = 1000;

% SSA data in folders gamma[i], i refers to gamma value
SSAfol = {'gamma5', 'gamma10', 'gamma20', 'gamma30'};
lenSSAfols = length(SSAfol);

% QB divided data in folders qb[i], i repeats, for both qbdiv and qbsave
QBfol = {'qb1', 'qb2', 'qb3'};
lenQBfols = length(QBfol);

% Root folder names
folder.IFS = 'processed ifs'; % save processed ifs data
folder.SSAload = 'ssa data/'; % ssa data set
folder.QBload = 'qb data/qb div/'; % divided qb data set
folder.QBsave = 'qb data/qb save/'; % save processed qb data

% File name roots for saving generated data or loading existing data
fileRoot.IFS = 'ifs';
fileRoot.loadQB = 'QBdiv';
fileRoot.saveQB = 'QBpcs';

% Temporary folder that will be reset each loop iteration
folTemp = folder;
thisDir = cd;


% Main folders to consider - sort into SSAfol names
cd(folder.SSAload);
repSSA = dir('ssa*');
nIndep = length(repSSA); % no. repeats
cd(thisDir);

for ind = 1:nIndep
    
    % Loop variables
    init = cell(lenSSAfols, lenQBfols);
    beta = init; gamma = init;
    mseEst = init; mseEst2 = init; meanEst = init;
    qbArea = init; numPreal = init; numPest = init;
    
    %% Main loops: runIntegrateFire over beta, repeat on qb, then for gamma
    
    for i = 1:lenSSAfols
                
        % For each gamma do this
        folTemp.SSAload = [folder.SSAload '/' repSSA(ind).name '/' SSAfol{i}];
        disp(['Processing SSA: ' folTemp.SSAload]);
        
        for j = 1:lenQBfols
            % For each qb repeated set do this
            folTemp.QBload = [folder.QBload QBfol{j}];
            folTemp.QBsave = [folder.QBsave QBfol{j}];
            disp(['Processing QB: ' folTemp.QBload]);
            
            % Main code that runs IFS on beta for a given gamma
            ifs = runIntegrateFire(folTemp, fileRoot, runs, runsdiv);
            
            % Parameters for simulations
            gamma{i, j} = ifs.gamma;
            beta{i, j} = ifs.beta;
            
            % MSE values
            mseEst2{i, j} = ifs.mseEst;
            mseEst{i, j} = ifs.mseEst;
            meanEst{i, j} = ifs.meanEst;
            
            % Integrate-fire data
            qbArea{i, j} = ifs.qbArea;
            numPest{i, j} = ifs.numPest;
            numPreal{i, j} = ifs.numPreal;
            
            % Order the beta and corresponding variables
            [beta{i, j}, sortID] = sort(beta{i, j});
            mseEst{i, j} = mseEst{i, j}(sortID);
            meanEst{i, j} = meanEst{i, j}(sortID);
            mseEst2{i, j} = mseEst2{i, j}(sortID);
            qbArea{i, j} = qbArea{i, j}(sortID);
            numPest{i, j} = numPest{i, j}(sortID);
            numPreal{i, j} = numPreal{i, j}(sortID);
            
            % Save temporary progress
            cd(folder.IFS);
            progName = [QBfol{j} '_' num2str(i) '.mat'];
            save(progName);
            cd(thisDir);
        end
        
        % Save progress on each gamma
        gammaName = [SSAfol{i} '.mat'];
        save(gammaName);
        disp('***************************************************************');
        disp(['Completed: ' num2str(i) ' of ' num2str(lenSSAfols)]);
        disp('***************************************************************');
    end
    
    % Run time
    tbatch = toc/60;
    disp(['Simulation run time = ' num2str(tbatch) ' mins']);
    
    % Get all m files called in creating this simulated data
    currName = dbstack; % gives a struct with current m file name
    currName = currName.file;
    [fList,pList] = matlab.codetools.requiredFilesAndProducts(currName);
    % Save all data with function list
    save(['complete_ ' num2str(ind) '_' num2str(lenSSAfols) '_' num2str(lenQBfols)]);
    
    %% Analysis and plotting for each gamma
    
    % Plot MSE across gammas on separate subplots
    figure;
    for i = 1:lenSSAfols
        subplot(ceil(lenSSAfols/2), 2, i);
        hold all
        for j = 1:lenQBfols
            plot(beta{i, j}, mseEst{i, j}, '-') %beta{i, j}, mseEst2{i, j}, 'o');
        end
        hold off
        xlabel('\beta');
        ylabel('MSE');
        title(['\gamma = ' num2str(gamma{i, 1}(1))]);
    end
    % Plot mean error across gammas on separate subplots
    figure;
    for i = 1:lenSSAfols
        subplot(ceil(lenSSAfols/2), 2, i);
        hold all
        for j = 1:lenQBfols
            plot(beta{i, j}, meanEst{i, j});
        end
        hold off
        xlabel('\beta');
        ylabel('mean error');
        title(['\gamma = ' num2str(gamma{i, 1}(1))]);
    end
    % QB area (average)
    figure;
    for i = 1:lenSSAfols
        subplot(ceil(lenSSAfols/2), 2, i);
        hold all
        for j = 1:lenQBfols
            plot(beta{i, j}, qbArea{i, j});
        end
        hold off
        xlabel('\beta');
        ylabel('avg QB area');
        title(['\gamma = ' num2str(gamma{i, 1}(1))]);
    end
end