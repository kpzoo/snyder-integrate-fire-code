% Function to load SSA data, combine with QBs and then run IFS
function outIFS = combineDataIFS(SSAname, runs, runsdiv, folder, fileRoot)

% Assumptions and modifications
% - based on threshPeakDetectFn3
% - getQBSSAStream assumes QBdiv[num] and gives QBpcs[num] and data folder
% - assumes process of interest is an IPP, so 2 state MC
% - works with single SSA data file
% - only calculates Snyder on the estimated photons
% - calls functions to estimate photons with integrate-fire

% Home directory
thisDir = cd;

% Controlling inputs
%loadname = 'QBreorg';
%savename = 'QBsaveSum';
%loadname = 'qb div';
%savename = 'qb save';

% QB stream inputs
nStep = 80;
%runsdiv = 1000;

%% Combine the SSA data with QBs by overlapping based on x2 event times

% Load Gillespie SSA data and extract variables
cd(folder.SSAload);
gillData = load(SSAname, 'outGil', 'params', 'Q');
cd(thisDir);

% Molecular count data 
X = gillData.outGil.X;
T = gillData.outGil.T;
x1 = X(:, 1); x2 = X(:, 2);

% Remove transient
T = T - T(1);
x2 = x2 - x2(1);

% Simulation parameters for SSA with gamma to control QB combination
beta = gillData.params.beta;
gamma = gillData.params.gamma;
alpha = gillData.params.alpha;
disp(['Loaded SSA: [beta gamma] = ' [num2str(beta) ' ' num2str(gamma)]]);

% Load variables from the Nikolic biological data and obtain an overlapped
% stream of QBs based on specified runs
outStream = getQBSSAStream(T, x2, x1, runs, runsdiv, nStep, folder, fileRoot);
Tset = outStream.Tset;
Iset = outStream.Iset;
%Tevent = outStream.Tevent;

% Interpolate the intensity function as did for QBs
[lamset , ~] = getSamplesExtrapQuad(alpha*x1, T, Tset, NaN);
% Remove time series entries for which there is no more intensity data
idrem = find(isnan(lamset));
idkeep = setdiff(1:length(Tset), idrem);
lamset = lamset(idkeep); % as don't want QBs beyond intensity range

% If removing these nan values then need to also remove the corresponding
% values in the overlapped stream (which was shifted with x2 events)
Tset = Tset(idkeep);
Iset = Iset(idkeep);

% Ensure that Tset has even intervals
dTset = diff(Tset);
if ~all(abs(dTset - dTset(1)) < 10^-8)
    assignin('base', 'dTset', dTset);
    error('Tset has uneven spacing');
end
clear('dTset');

%% Convert overlapped stream into estimated photons and filter

% Get all event times of true photons
Pset = outStream.Tevent;

% Get the estimated integrate fire photons at Test
[Test, dist, thresh] = getIntegFireResults(Tset, Iset, Pset, 0);
Test = Test';
clear('Tset', 'Iset');
% Counts of real and estimated photons
lenTreal = length(Pset);
lenTest = length(Test);

% State limit information
Slim = gillData.params.SlimSet;
S = diag(Slim.min(1):Slim.max(1));
lenS = length(S);

% Uniform (discrete) prior and Lam for Snyder solution
q0 = (1/lenS)*ones(1, lenS);
Lam = alpha*S;

% Regenerate datasets to remove repetitions and make a new [T X] set
%[TestNew , ~] = reconstructTX(T, [x1 x2], Test);
%[TrealNew , ~] = reconstructTX(T, [x1 x2], Pset);

% Use restricted SSA data and generate new form with estimated photons
%[TestNew, XestNew] = reconstructTX(outStream.Tres, [outStream.x1res outStream.x2res], Test);
%[TNew, XNew] = reconstructTX(outStream.Tres, [outStream.x1res outStream.x2res], Pset);

% Main Snyder filter code
outFil = filterDSPPIden(Test, q0, gillData.Q, Lam, lenS);

% Check the probability vectors at event times
if max(abs(sum(outFil.qev, 2) - 1)) > 10^-9
    error('The probability density sum is not close enough to 1');
end
if min(min(outFil.qev)) < 0
    error('The probability density has negative values');
end

% Posteriors from Snyder filter across time
Tfil = outFil.Tset;
Qfil = outFil.Qset;
Qfilgood = cell(1, 1);
Tfilgood = cell(1, 1);

% Remove indices with 1 element or less (from redundancy)
innerLen = zeros(1, length(Qfil));
j = 1;
for i = 1:length(Tfil)
    innerLen(i) = length(Tfil{i});
    if innerLen(i) > 1
        % Assign elements with > 1 length only
        Qfilgood{j} = Qfil{i};
        Tfilgood{j} = Tfil{i};
        j = j + 1;
    end
end
% Lengths <= 1 cannot be interpolated and useless for MSE
goodid = find(innerLen > 1);
lenGood = length(goodid);
disp(['No. bad values = ' num2str(length(Tfil) - lenGood)]);

% Reassign names and remove excess variables
Qfil = Qfilgood;
Tfil = Tfilgood;
clear('Qfilgood', 'Tfilgood');

% Get corresponding x1 estimates
x1fil = cell(1, 1);
j = 1;
for i = 1:length(Qfil)
    % x1 estimates are posteriors times state space
    x1fil{j} = Qfil{i}*diag(S);
    j = j + 1;
end
    
% Get MSE and clear extra variables
x1stats = getStatsMSE2(Tfil, x1fil, T, x1);

% Convert all cells to a single vectors
tcat = Tfil{1};
qcat = Qfil{1}; 
xcat = x1fil{1};
for i = 2:length(Tfil)
    % Concatenate the arrays in each cell
    tcat = cat(1, tcat, Tfil{i});
    qcat = cat(1, qcat, Qfil{i});
    xcat = cat(1, xcat, x1fil{i});
end
x1stats2 = getStatsAltMSE(tcat, xcat, x1, T);

cellvars = {'Qfil', 'Tfil', 'x1fil', 'tcat', 'xcat', 'qcat'};
clear(cellvars{:});

% Save simulation parameters
outIFS.beta = beta;
outIFS.gamma = gamma;
outIFS.avgQBArea = thresh;
outIFS.dist = dist;

% Estimated and true photons
outIFS.Test = Test;
outIFS.Treal = Pset;
outIFS.lenTreal = lenTreal;
outIFS.lenTest = lenTest;
%outIFS.estGill = [TestNew, XestNew];
%outIFS.resGill = [TNew, XNew];

% MSE values
outIFS.meanEst = x1stats.vals(1);
outIFS.mseEst = x1stats.vals(3);
outIFS.mseEst2 = x1stats2.vals(3);
