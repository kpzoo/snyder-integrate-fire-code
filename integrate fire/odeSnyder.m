% Function to solve ODEs on the Snyder filter between x2 birth events -
% uses Q matrix and the hybrid stochastic method but keep dq/dt vs dq/dI
function dy = odeSnyder(ts, y, Q, Lam)

% Assumptions and modifications
% - y is a column vector
% - diagS is the diagonal of the state space S matrix

% Lamhat, the estimated Lam constant
Lamhat = y'*diag(Lam);
Lamhat = Lamhat*eye(length(Lam));

% Snyder ODE set for Markov jump process
dy = y'*(Q - Lam + Lamhat);
% Output must be column vector
dy = dy';


