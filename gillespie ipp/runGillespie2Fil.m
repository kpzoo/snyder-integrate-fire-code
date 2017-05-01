% Function to run Gillespie simulations and Snyder filtering as required on
% multi-state MCs with complicated stationary distributions
function params = runGillespie2Fil(alpha, k, simname, stateMax, jumpMax, skip_Sny, folStore)

% Assumptions and Modifications
% - skips SSA, does filtering only
% - added path for saving
% - added control of save folders and skip_Sny input
% - handles general Markov chains with non-tridiagonal intensity matrices
% - Gillespie SSA done via dsppMC code
% - assumes proportional relation x2dot = alpha*x1
% - loops across beta for a fixed state transition r_const term
% - inserts zeros for missed reaction jump values e.g. if 1 and 3 jumps
% given rates then puts a 0 rate for jumps of size 2
% - can skip_SSA to allow for simply Snyder simulations
% - bulk reactions mean non-nearest neighbour transitions


% Set boolean to control the saving of data and figures 
save_sim = 1;

% Mode to skip SSA or Snyder simulations if want to reuse data
skip_SSA = 1;
skip_Snyder = skip_Sny;

%% Gillespie simulation to obtain counts in time

% This cell is solely concerned with the observed process and provides the
% actual modulated intensity that is to be estimated

% Run SSA if boolean set else run filter with pre-existing simulation data
if ~skip_SSA
    % Assign state space limits as [MCstate obs counts]
    params.SlimSet.max = [stateMax inf];
    % Min states assumed as 0
    params.SlimSet.min = [0 0];

    % Set appropriate r_const for specified jumpMax, first 2 elements are
    % the nearest neighbour transitions (r_const is rate constants)
    switch(jumpMax)
        case 1
            % IPP model
            r_const = [k k alpha];
            inpQ = [r_const(2) 0 0];
        case 16
            % Bimodal 16 state
            r_const = [k k zeros(1, 2*jumpMax-4) k k alpha];
            inpQ = [r_const(2) r_const(2*jumpMax) 0]; 
        otherwise
            error('Value of jumpMax not supported');
    end

    % Length of r_const should be twice jumpMax and a single x2 reaction
    if length(r_const) ~= 2*jumpMax + 1
        error('Inconsistent r_const definition for specified jumpMax');
    end

    % Here transit represents the incremental changes that 
    % bulk reactions yield - odd reactions are births and even deaths
    [molecType, crossType, reacType, bulk, transit] = setMCParams(jumpMax, r_const, 1);
    % Assign outputs to structure
    params.molecType = molecType; params.crossType = crossType;
    params.reacType = reacType; params.bulk = bulk;
    params.transit = transit; params.r_const = r_const;

    % Check that bulk matches the jumpMax setting and lengths are correct
    if max(params.bulk) ~= jumpMax
        error('Inconsistency between bulk and jumpMax');
    end
    % Ensure all lengths are consistent
    if all([length(params.bulk) length(params.r_const) length(params.transit)...
            length(params.molecType) length(params.crossType)] ~= length(params.reacType))
        error('The lengths of the input simulation vectors do not match the number of reactions');
    else
        nReacs = length(params.reacType);
    end

    % Assign simulation control parameters and initial populations
    params.Nstart = 5000;
    params.N = 30000;
    params.len = 2; % 2 entries as x1 and x2 (observed counts)
    params.x0 = zeros(1, params.len); 
    params.nReacs = nReacs;

    % Assuming proportional structure between x2dot and x1 
    params.kgain = params.r_const(end);
    params.coeff = [params.kgain 0];

    % Assumes 1st state is 0 and gives bimodal or IPP Q matrix
    if params.SlimSet.min(1) == 0
        % No. of states in space
        nState = length(params.SlimSet.min(1):params.SlimSet.max(1));
        % Bimodal vs uniform distribution on multi-state models
        statdistr = 2;
        % Calculate Q matrix and stationary distribution Pi
        [Q, P, Pi] = setQMx(nState, statdistr, inpQ);
    else
        error('Q matrix code does not support non-zero minimum states');
    end

    % Parameters for Gillespie SSA - Q matrix and intensity settings
    params.Q = Q; params.P = P; params.Pi = Pi;
    params.alpha = alpha;
    params.gamma = 1/(100*k);
    params.beta = alpha/k;

    % Run the Gillespie simulations and extract the actual rate process and
    % counting observations (in X and T)
    disp('Simulation started');
    outGil = dsppMCQ(params);
    disp('Simulation complete');

    % Check the rates match those delimited as possible by Q
    checkQ = checkRateQ(params.Q, outGil.Xdot);
    if ~all(checkQ == 1)
        error('Some rates do not match Q values');
    end
    % Check unique X values match those of the space in number
    X = outGil.X;
    actualStates = unique(X(:, 1));
    if length(actualStates) ~= length(params.SlimSet.min(1):params.SlimSet.max(1))
        disp('The simulation has not traversed all states within Nstart to N');
    end

    % Make a new directory and save Gillespie data
    thisDir = cd;
    cd('gil data');
    if exist(folStore, 'dir') ~= 7
        % Folder does not exist unless 7 returned
        mkdir(folStore);
    end
    cd(folStore);
    disp(['Saved file: ' simname]);
    save(simname);
    cd(thisDir);
end

%% Snyder filter the data at the front end of the cascade (no biology)

% Run Snyder filter only if boolean set
if ~skip_Snyder
    % Load data file and clear window
    clc
    load(simname);

    % Obtain time and remove offset from the transient samples
    T = outGil.T;
    %Tref = T(1);
    %T = T - Tref;

    % Obtain molecular counts
    lenT = length(T);
    X = outGil.X;
    x1 = X(1:lenT, 1);
    x2 = X(1:lenT, 2);

    % Obtain the state space of x1 and check dimensions
    SlimSet = params.SlimSet;
    S = diag(SlimSet.min(1):SlimSet.max(1));
    lenS = length(S);
    
    
    % Uniform (discrete) prior and Lam for Snyder solution
    q0 = (1/lenS)*ones(1, lenS);
    Lam = params.alpha*S;
    
    % Get event times of observed process x2
    dx2 = [0; diff(x2)]; % 0 to get correct index
    % Indices of points where x2 increments
    Tevent = T(dx2 == 1);

    % Apply suitable inputs to the filter and extract outputs
    disp('Running Snyder filter');
    outFil = filterDSPP(Tevent, q0, params.Q, Lam, lenS);
    disp('Finished estimate of modulating intensity');

    % Check the probability vectors at event times
    if max(abs(sum(outFil.qev, 2) - 1)) > 10^-9
        error('The probability density sum is not close enough to 1');
    end
    if min(min(outFil.qev)) < 0
        error('The probability density has negative values');
    end

    % Obtain parameters for plots and statistics
    Tset = outFil.Tset;
    Qset = outFil.Qset;
    x1Set = cell(1, 1);
    for i = 2:length(Qset)
        % x1 estimates are posteriors times state space
        x1Set{i} = Qset{i}*diag(S);
    end
    
    % Convert all cells to a single vectors
    tcat = Tset{2};
    qcat = Qset{2}; % first element is always empty
    xcat = x1Set{2};
    for i = 3:length(Qset)
        % Concatenate the arrays in each cell 
        tcat = cat(1, tcat, Tset{i}); 
        qcat = cat(1, qcat, Qset{i});
        xcat = cat(1, xcat, x1Set{i});
    end
     
    % Calculate MSE from posteriors
    params.x1stats1 = getStatsMSE(Tset, x1Set, T, x1, Tevent);
    disp(params.x1stats1);
    params.x1stats2 = getStatsAltMSE(tcat, xcat, x1, T);
    disp(params.x1stats2);
    
    % Clear variables not to be saved
    cellvars = {'Qset', 'Tset', 'x1Set', 'tcat', 'xcat', 'qcat'};
    clear(cellvars{:});

    % Save filtered data if boolean set
    if save_sim
        cd('fil data');
        if exist(folStore, 'dir') ~= 7
            % Folder does not exist unless 7 returned
            mkdir(folStore);
        end
        cd(folStore);
        filname = ['filter_' simname];
        save(filname);
        cd(thisDir);
    end
    disp('Simulation and post processing complete');

end