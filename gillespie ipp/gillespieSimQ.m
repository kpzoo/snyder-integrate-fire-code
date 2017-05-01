% Stochastic simulation algorithm of Gillespie for a process with 2
% molecular species, x1 (hidden MC) and x2 (observed DSPP)
function [X, Xdot, T] = gillespieSimQ(inpGill)

% Assumptions and modifications
% - x2 assumed to only have births with rate lam = alpha*x1
% - x1 modulates the intensity of x2, no other cross-relationships
% - allows non-nearest neighbour reactions
% - transit matrix corresponds to r_const form [bir dea bir dea bir dea...]
% - reac rates calculated in this order, deaType are even, birType odd idx
% - used a lookup on the Q matrix

% Initialisation of SSA, start and total event count
len = inpGill.len; x0 = inpGill.x0;
N = inpGill.N; Nstart = inpGill.Nstart;
Q = inpGill.Q;

% Rate constants and transitions
r_const = inpGill.r_const;
transit = inpGill.transit;

% Set initial parameters for time and populations
t = zeros(N, 1);
z = zeros(N, len);
z(1, :) = x0;

% Set initial rates and reaction incrementation - based on r_const
nReacs = length(r_const);
alpha = zeros(N, nReacs);

% Loop across specified number of iterations and simulate
for i = 2:N
    % Update molecular numbers and other data sequentially
    x = z(i - 1, :);
    told = t(i-1);
    
    % Lookup x1 birth and death rates across its chain from Q
    rdotx1 = getx1Rates(Q, x(1), nReacs-1);
    % Rate of x2 births is lam = alpha*x1, alpha should be r_const(end) 
    rdotx2 = r_const(end)*x(1);
    rdot = [rdotx1 rdotx2];

    % Get exponential time to next reaction
    rdotsum = sum(rdot);
    tnex = told - log(rand)/rdotsum;
    % Find the next reaction based on propensities (rate ratios)
    rdot_ratio = rdot/rdotsum;
    reac = 1 + sum(rand > cumsum(rdot_ratio)); % choice via uniform distr
    % Append molecular count based on reaction of choice
    try
        xnex = x + transit(:, reac)';
    catch
        assignin('base', 'reac', reac);
        return;
    end
    % Assign variables
    try
        alpha(i-1, :) = rdot;
    catch
        assignin('base', 'rdot', rdot);
        assignin('base', 'alpha', alpha(i-1, :));
        error('Dimension mismatch');
    end
    t(i) = tnex;
    z(i, :) = xnex;

    % Catch a possible error
    if rdotsum == 0
        assignin('base', 'z', z(1:i, :));
        assignin('base', 'alpha', alpha(1:i-1, :));
        error(['All rates are zero at i = ' num2str(i)]);
    end

end

% Save data from simulation for post processing and account for the control
x = z;
xdot = alpha;

% Account for equilibrium and remove last value (in previous code would
% calculate the next step at i = N)
X = x(Nstart:N-1, :);
Xdot = xdot(Nstart:N-1, :);
T = t(Nstart:N-1, :);
