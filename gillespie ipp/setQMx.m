% Function to construct a Q matrix for a given stationary state distribution
% Pi, the code assumes that Pi has no zero entries
function [Q, P, Pi] = setQMx(nState, statdistr, inpQ)

% Assumptions and modifications
% - assumes symmetric IPP and bimodal models only
% - inpQ is an input that controls the Q formulation
% - assumes lowest state is 0
% - constructs a tridiagonal Q based on vector of birth and death rates
% - adds extra entries for modal shifts

if nState == 2
    % IPP model so set Q matrix directly
    if length(inpQ(inpQ > 0)) > 1
        % Only 1 value must be non-zero, inpQ(1)
        error('Incorrect input, just need symmetric gain k');
    else
        % Simple Q matrix for symmetric rates
        Q = inpQ(1)*[-1 1; 1 -1];
        Pi = [0.5 0.5];
        disp('Symmetric IPP model specified');
    end  
else
    % State space for model with > 2 states
    states = 0:nState-1;
    
    % Construct different desired stationary distributions
    switch(statdistr)
        case 1
            % Discrete uniform distribution
            Pi = ones(1, nState)/nState;
        case 2
            % Bimodal distribution obtained from 2 Gaussian distributions
            if nState <= 5
                error('Not enough states for sensible bimodal form');
            end
            % Location and spread of Gaussians
            mu1 = nState/4; sigma1 = mu1/4;
            mu2 = 3*nState/4; sigma2 = mu2/12;
            % Bimodal stationary distribution (discrete)
            Pi = normpdf(states, mu1, sigma1) + normpdf(states, mu2, sigma2);
        otherwise
            error(['No distribution available for statdistr = ' num2str(statdistr)]);
    end
    
    % Ensure Pi sums to 1
    Pi = Pi/sum(Pi);
    % Ensure no zero or negative Pi entries
    if any(Pi <= 0)
        error('No support for Pi with non-positive entries');
    end
    
    % Test input size and assign the single gain for deaths in
    % addition to 1 value for eps's (mode switch) gain
    nExp = 3;
    if length(inpQ) ~= nExp
        disp(['Expected ' num2str(nExp) ' inputs']);
        error('Incorrect unconstrained input size for case 5');
    else
        % Death rate vector
        b = inpQ(1)*states(2:end);
        % Modal switch rate
        eps1 = inpQ(end-1)*states(end);
        % Difference between switch rates, must be 0 for this case
        delta = inpQ(end);
    end
    
    % Matrix with single entries to allow transition between the modes
    % of the bimodal setup - assuming 2 max values 
    maxStates = states(Pi == max(Pi));
    maxid = maxStates + 1;
    if length(maxid) ~= 2
        error('Bimodal distribution does not give 2 maxima');
    end
    
    % Calculate second epsilon value
    %eps2 = (delta + eps1*Pi(maxid(1)))/Pi(maxid(2));
    eps2 = eps1;
    
    % By allowing delta = 0 a homogeneous solution is obtained for the
    % birth rate vector a from the specified b and eps's
    if delta == 0
        r = Pi(1:end-1)./Pi(2:end);
        % Birth rate vector
        a = b./r;
        % Construct the bimodal Q
        Q = diag(a, 1) + diag(b, -1);
        Q(maxid(1), maxid(2)) = eps1;
        Q(maxid(2), maxid(1)) = eps2;
        Q = Q - diag(sum(Q, 2));
    else
        error('Non zero delta methods not supported');
    end 
    disp([num2str(nState) ' state model specified']);
end

% Check that the Q matrix is properly conditioned
if any(diag(Q) > 0)
    error('Q matrix has a positive diagonal element');
end
sumRow = sum(sum(Q, 2));
if sumRow > 10^-10
    error('The row sums are not close enough to 0');
end
Qnodiag = Q - diag(diag(Q));
if min(min(Qnodiag)) < 0
    error('Off diagonal elements are not all non-negative');
end
% Check that the Q really has Pi as its null space
nu = Pi*Q;
sum_nu = sum(nu);
if sum_nu < 10^-10
    disp('Estimated Q matrix seems correct');
end
% Calculate P matrix from Kolmogorov solution
t = 100000;
P = expm(t*Q);

% Check that time is high enough
eP = abs(P^2 - P);
if max(max(eP)) > 10^-8
    warning('Mat:Pcalc', 'Specified time not large enough in P calculation');
end