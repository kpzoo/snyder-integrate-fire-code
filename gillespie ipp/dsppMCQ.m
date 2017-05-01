% Function to generate a general Markov Chain Gillespie simulation
function output = dsppMCQ(params)

% Assumptions and modification
% - works with a specified Q matrix
% - performs SSA

% Obtain input parameters - no. molecules, iterations, rate extrema and
% molecular ICs, start of equilibrium in terms of N
len = params.len; x0 = params.x0;
N = params.N; Nstart = params.Nstart;

% Rate structure based inputs - odd reactions are births and even deaths,
% bulk is integer change of reactionss
molecType = params.molecType;
crossType = params.crossType;

% Rate constants and transition matrix for changes in molecular numbers
transit = params.transit;
r_const = params.r_const;

% Q matrix assignments
Q = params.Q;
inpGill.Q = Q;

% Inputs to Gillespie Markov algorithm
inpGill.len = len; inpGill.x0 = x0;
inpGill.N = N; inpGill.Nstart = Nstart;
inpGill.r_const = r_const; inpGill.transit = transit; 

% Ensure molecular types sensible
if ~all(ismember([molecType crossType], 1:len))
    assignin('base', 'molecType', molecType);
    assignin('base', 'len', len);
    error('Incorrect molecular identifiers specified');
end

% Run the Gillespie algorithm for specified Markov chain
disp('Gillespie simulation of Markov chain started');
[X, Xdot, T] = gillespieSimQ(inpGill);
disp('Gillespie simulation of Markov chain complete');
% Assign outputs to a single data structure
output.X = X; output.Xdot = Xdot; output.T = T;
