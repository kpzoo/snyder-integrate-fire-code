% Function implements filtering equations from Snyder in combination with
% use of ODE integration
function outFil = filterDSPP(Tevent, q0, Q, Lam, lenS)

% Assumptions and modifications
% - x2 is the observed process (molecular counts)
% - T is the time vector for x1 and x2 reactions of Gillespie SSA form

% Length of event stream
nEv = length(Tevent);

% Declare variables and initialise
qev = zeros(nEv, lenS);
qev(1, :) = q0;
tev = zeros(nEv, 1);

% Cell to save output of ODE solver and set options
Qset = cell(1, 1);
Tset = cell(1, 1);
options = odeset('NonNegative', 1:lenS);
% options = odeset(options, 'Refine', 20); % <----- this option works well
% options = odeset(options, 'RelTol', 0.001);

% Loop across events integrating the Snyder equations with the ODE solver
% and then applying jump corrections due to x2 events
for i = 2:nEv
    % Solve ODEs continuously with setting of options, integration limits
    [tset, qset] = ode113(@(ts, y) odeSnyder(ts, y, Q, Lam),...
        [Tevent(i-1) Tevent(i)], qev(i-1, :)', options);

    % Assign output value at event times
    qev(i, :) = qset(end, :);
    tev(i) = tset(end);
    % Distributions across all ODE evaluation times
    Qset{i} = qset;
    Tset{i} = tset;

    % Check raw ODE q output
    if any(qev(i, :) < -10^-8)
        assignin('base', 'qODE', qev(i, :));
        error(['qODE distribution has negative entries at i =' num2str(i)]);
    end
    if max(abs(sum(qev(i, :)) - 1)) > 10^-4
        assignin('base', 'qODE', qev(i, :));
        error(['qODE distribution does not sum to 1 at i = ' num2str(i)]);
    end

    % Calculate perturbation on q due to jump
    pert = diag(Lam)';
    qev(i, :) = qev(i, :).*pert./(qev(i, :)*pert');

    % Display progress and debug set of results from ODE solver
    disp(['Finished iteration: ' num2str(i-1) ' of ' num2str(nEv)]);

end

% Assign output data to a single structure
outFil.qev = qev;
outFil.t = tev;
outFil.Tevent = Tevent;
outFil.Qset = Qset;
outFil.Tset = Tset;
