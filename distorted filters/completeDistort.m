% Run deterministic delayed photons on repeat streams at each gamma
clearvars
clc
close all

% Assumptions and modifications
% - assume the distorted photons are already available
% - assume a folder structure for repeat photon sets
% - loop 1 is on gamma, 2 is on photons and 3 on beta

%% Setup the problem and folder structure

tic; 

% Options for deterministic delay
tau = 43.3;
limPhoton = 8000; % photon limit

% SSA data in folders gamma[i], i refers to gamma value
SSAfol = {'gamma5', 'gamma10', 'gamma20', 'gamma30'};
lenSSAfols = length(SSAfol);

% Root folder names
folder.det = 'processed det'; % save processed delayed data
folder.SSAload = 'ssa data/'; % ssa data set

% File name roots for saving generated data or loading existing data
fileRoot.det = 'det';

% Temporary folder that will be reset each loop iteration
folTemp = folder;
thisDir = cd;

% Get Gillespise SSA source file length
cd([folder.SSAload SSAfol{1}]);
gillFiles = dir('*.mat');
cd(thisDir);
% No. files to process (assumed same in each folder)
flen = length(gillFiles); 

% Loop variables
init = zeros(lenSSAfols, flen);
beta = init; gamma = init; alpha = init;
mseEst = init; mseEst2 = init; meanEst = init;

%% Main loops: runIntegrateFire over beta, repeat on qb, then for gamma

for j = 1:lenSSAfols
   % For each gamma do this 
    folTemp.SSAload = [folder.SSAload SSAfol{j}];
    disp(['Processing SSA: ' folTemp.SSAload]);
    
    % List of data in each gamma folder
    cd([folder.SSAload SSAfol{1}]);
    gillFiles = dir('*.mat');
    cd(thisDir);
    if length(gillFiles) ~= flen
        error('No. files not identical for each gamma folder');
    end
    
    % Load each SSA file and process/filter
    for i = 1:flen
        % SSA data to be loaded
        cd([folder.SSAload SSAfol{j}]);
        gillFiles = dir('*.mat');
        gillDataName = gillFiles(i).name;
        gillData = load(gillDataName, 'outGil', 'params', 'Q');
        cd(thisDir);
        
        % Extract SSA molecular counts and event times and state space
        X = gillData.outGil.X;
        T = gillData.outGil.T;
        Slim = gillData.params.SlimSet;
        Q = gillData.Q;
        
        % Get 3 key simulation settings
        beta(j, i) = gillData.params.beta;
        gamma(j, i) = gillData.params.gamma;
        alpha(j, i) = gillData.params.alpha;
        
        % Remove rest of loaded data
        clearvars('gillData');
        
        % Prepare the SSA data by removing offsets on x2 and T due to transient
        X(:, 2) = X(:, 2) - X(1, 2);
        T = T - T(1);
        
        % Limit the photon stream based on limPhoton input
        if X(end, 2) < limPhoton
            warning('Mat:fewPh', ['There are less than ' num2str(limPhoton) ' photons present ']);
        else
            disp(['Clipping Gillespie data to ' num2str(limPhoton) ' photons']);
            idlim = find(X(:, 2) == limPhoton, 1, 'last');
            X = X(1:idlim, :);
            T = T(1:idlim);
        end
        
        % Get event times of observed process (photon times)
        dx2 = [0; diff(X(:, 2))]; % 0 to get correct index
        % Indices of points where new photon comes
        Tevent = T(dx2 == 1);
        
        % Get state space and uniform prior
        S = diag(Slim.min(1):Slim.max(1));
        lenS = length(S);
        
        % Uniform (discrete) prior and Lam for Snyder solution
        q0 = (1/lenS)*ones(1, lenS);
        Lam = alpha(j, i)*S;
        
        % Add deterministic delay and use optimal delayed Snyder filter
        outDet = detDelayProcessing(Tevent, tau, T, X, q0, Q, Lam, lenS, S);
        
        % Get MSE
        meanEst(j, i) = outDet.meanEst;
        mseEst(j, i) = outDet.mseEst;
        mseEst2(j, i) = outDet.mseEst2;
    end

        
    % Save progress on each gamma
    gammaName = [SSAfol{j} '.mat'];
    save(gammaName);
    disp('***************************************************************');
    disp(['Completed: ' num2str(j) ' of ' num2str(lenSSAfols)]);
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
cd(folder.det);
save(['completeDet_ ' num2str(lenSSAfols)]);
cd(thisDir);

%% Analysis and plotting for each gamma

% Sort all the data
for i = 1:lenSSAfols
    [beta(i, :), sortID] = sort(beta(i, :));
    mseEst(i, :) = mseEst(i, sortID);
    meanEst(i, :) = meanEst(i, sortID);
    mseEst2(i, :) = mseEst2(i, sortID);
end

% Plot MSE across gammas on separate subplots
figure;
for i = 1:lenSSAfols
    subplot(ceil(lenSSAfols/2), 2, i);
    hold all
    plot(beta', mseEst', '-', beta', mseEst2', 'o');
    hold off
    xlabel('\beta');
    ylabel('MSE');
    title(['\gamma = ' num2str(gamma(i, 1))]);
end
% Plot mean error across gammas on separate subplots
figure;
for i = 1:lenSSAfols
    subplot(ceil(lenSSAfols/2), 2, i);
    hold all
    plot(beta', meanEst');
    hold off
    xlabel('\beta');
    ylabel('mean error');
    title(['\gamma = ' num2str(gamma(i, 1))]);
end

