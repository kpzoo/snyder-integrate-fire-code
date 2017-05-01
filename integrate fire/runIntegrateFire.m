% Function that runs integrate-fire-Snyder on all the files in a given
% subfolder and returns estimated photons and MSE values
function ifs = runIntegrateFire(folder, fileRoot, runs, runsdiv)

% Assumptions and modifications
% - a single threshold given by the average QB area (local) is used

% Existing directory
thisDir = cd;

% Obtain source files - Gillespie simulation data
cd(folder.SSAload);
files = dir('*.mat');
cd(thisDir);
% No. files to process
flen = length(files);

% Simulation variables
betaSet = zeros(1, flen);
gammaSet = zeros(1, flen);
avgQBArea = zeros(1, flen);

% MSE and metric variables
meanEst = zeros(1, flen);
mseEst = zeros(1, flen);
mseEst2 = zeros(1, flen);

% No. estimated and real photons in file
numPest = zeros(1, flen);
numPreal = zeros(1, flen);

% Loop across batch runs and obtain parameters of interest
for i = 1:flen
    % File to be processed by IFS
    SSAname = files(i).name;
    % Main function to combine SSA-QB data and run Snyder
    outIFS = combineDataIFS(SSAname, runs, runsdiv, folder, fileRoot);
    
    % Simulation parameters and threshold of avg QB area
    betaSet(i) = outIFS.beta;
    gammaSet(i) = outIFS.gamma;
    avgQBArea(i) = outIFS.avgQBArea;
    
    % Snyder MSE on estimated photons
    meanEst(i) = outIFS.meanEst;
    mseEst(i) = outIFS.mseEst;
    mseEst2(i) = outIFS.mseEst2;
    
    % Further variables to store how many photons are actually in the
    % stream and how much are produced by integrate-fire
    numPest(i) = outIFS.lenTest;
    numPreal(i) = outIFS.lenTreal;
    
    % Save and display progress
    % save(['IFSProg' num2str(i)]);
    disp('-------------------------------------------------------------');
    disp(['Completed iteration ' num2str(i) ' out of ' num2str(flen)]);
    disp('-------------------------------------------------------------');
end

% Data structure to output main results
ifs.beta = betaSet;
ifs.gamma = gammaSet;
ifs.qbArea = avgQBArea;
ifs.mseEst = mseEst;
ifs.mseEst2 = mseEst2;
ifs.meanEst = meanEst;
ifs.numPest = numPest;
ifs.numPreal = numPreal;

% Save completed data
cd(folder.IFS);
save(fileRoot.IFS); 
cd(thisDir);