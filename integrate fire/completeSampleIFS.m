% Run IFS code on several repeat photon streams at each gamma
clearvars
clc
close all

% Assumptions and modifications
% - expects 1 independent qb run per ssa run
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

% QB divided data in folders qb[i], i repeats, for both qbdiv and qbsave
cd(folder.QBload);
QBfol = dir('qb*');
cd(thisDir);
lenQBfols = 1; % because using 1 per ssa

% Main folders to consider - sort into SSAfol names
cd(folder.SSAload);
repSSA = dir('ssa*');
nIndep = length(repSSA); % no. repeats
cd(thisDir);

% Check have enough data
numQBfolsTot = length(QBfol);
if numQBfolsTot ~= nIndep
    error('Data files not of correct number');
end

% Create equivalent directories in qb save as qb div
cd(folder.QBsave);
for i = 1:numQBfolsTot
    % Create each folder if not present
    if exist(QBfol(i).name, 'dir') ~= 7
        % Folder does not exist unless 7 returned
        mkdir(QBfol(i).name);
    end
end
cd(thisDir);

% Time the code
tbatch = zeros(1, nIndep);


for ind = 1:nIndep
    
    % Loop variables
    init = cell(lenSSAfols, lenQBfols);
    beta = init; gamma = init;
    mseEst = init; mseEst2 = init; meanEst = init;
    qbArea = init; numPreal = init; numPest = init;
    
    %% Main loops: runIntegrateFire over beta, then for gamma, then repeat
    
    for i = 1:lenSSAfols
                
        % For each gamma do this
        folTemp.SSAload = [folder.SSAload repSSA(ind).name '/' SSAfol{i}];
        disp(['Processing SSA: ' folTemp.SSAload]);
        
        % For each qb repeated set do this
        folTemp.QBload = [folder.QBload QBfol(ind).name];
        folTemp.QBsave = [folder.QBsave QBfol(ind).name];
        disp(['Processing QB: ' folTemp.QBload]);
        
        % Main code that runs IFS on beta for a given gamma
        ifs = runIntegrateFire(folTemp, fileRoot, runs, runsdiv);
        
        % Parameters for simulations
        gamma{i} = ifs.gamma;
        beta{i} = ifs.beta;
        
        % MSE values
        mseEst2{i} = ifs.mseEst;
        mseEst{i} = ifs.mseEst;
        meanEst{i} = ifs.meanEst;
        
        % Integrate-fire data
        qbArea{i} = ifs.qbArea;
        numPest{i} = ifs.numPest;
        numPreal{i} = ifs.numPreal;
        
        % Order the beta and corresponding variables
        [beta{i}, sortID] = sort(beta{i});
        mseEst{i} = mseEst{i}(sortID);
        meanEst{i} = meanEst{i}(sortID);
        mseEst2{i} = mseEst2{i}(sortID);
        qbArea{i} = qbArea{i}(sortID);
        numPest{i} = numPest{i}(sortID);
        numPreal{i} = numPreal{i}(sortID);
        
        % Save temporary progress
        cd(folder.IFS);
        progName = [QBfol(ind).name '_' num2str(i) '.mat'];
        save(progName);
        cd(thisDir);
        
        % Save progress on each gamma
        disp('***************************************************************');
        disp(['Completed: ' num2str(i) ' of ' num2str(lenSSAfols)]);
        disp('***************************************************************');
    end
    
    % Run time
    if ind > 1
        tbatch(ind) = toc/60 - tbatch(ind-1);
    else
        tbatch(ind) = toc/60;
    end
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