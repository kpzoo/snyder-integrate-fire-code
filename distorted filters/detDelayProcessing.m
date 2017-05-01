% Function to deterministically delay and then optimally filter photons
function outDet = detDelayProcessing(Tevent, tau, T, X, q0, Q, Lam, lenS, S)

% Assumptions and modifications 
% - event times and delay (tau) in same units

% Deterministically delay the observed true phtoton times
Tdist = Tevent + tau;
outDet.tau = tau;

% Reconstruct an SSA type output with the distorted photons
[Tdel, Xdel] = getDistortedPhotonStream(T, X, Tdist);
outDet.Tdel = Tdel;
outDet.Xdel = Xdel;

% Apply standard Snyder filter to delayed data
outFil = filterDSPP(sort(Tdist), q0, Q, Lam, lenS);
qev = outFil.qev;

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

% Apply Kolmogorov exponential correction for posterior evolution of
% deterministically delayed events
expSol = expm(tau*Q);
for i = 1:length(Qfil)
    % Correction to optimally account for the delay
    Qfil{i} = Qfil{i}*expSol;
end
qev = qev*expSol;

% Get corresponding x1 estimates
x1fil = cell(1, 1);
j = 1;
for i = 1:length(Qfil)
    % x1 estimates are posteriors times state space
    x1fil{j} = Qfil{i}*diag(S);
    j = j + 1;
end
    
% Get MSE and clear extra variables
x1stats = getStatsMSE2(Tfil, x1fil, Tdel, Xdel(:, 1));

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
x1stats2 = getStatsAltMSE(tcat, xcat, Xdel(:, 1), Tdel);
cellvars = {'Qfil', 'Tfil', 'x1fil', 'tcat', 'xcat', 'qcat'};
clear(cellvars{:});


% Assign output structure
outDet.qev = qev;
outDet.Tdel = Tdel;
%outDet.Qfil = Qfil;
%outDet.Tfil = Tfil;
outDet.meanEst = x1stats.vals(1);
outDet.mseEst = x1stats.vals(3);
outDet.mseEst2 = x1stats2.vals(3);