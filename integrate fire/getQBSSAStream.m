% Optimised to reduce the number of saves and loads and to make use of more
% efficient interpolation as well as use nakeinterp1 and alters the
% savename and loadname style via assumption that QBs already divided

% Altered version of obtainStreamQB2 that works with several QB data sets
% in folder QBdata and combines them to obtain sufficient photons based on
% an input requirement - further no overlap functionality here

% Function to construct a QB stream based on the biological model's output
% and the x2 births with the option of overlapping QBs
function outComb = getQBSSAStream(T, x2, x1, runs, runsdiv, nStep, folder, fileRoot)

% Assumptions and modifications
% - from getQBstream3
% - assumes QB data already reorganised into runsdiv size files
% - ^should only reorganise once ever and fix? <--------------
% - QB data in folder structure: data -> loadname, name QBdiv[num] 
% - assumes all data at same sample frequency and size

% Set home directory
thisDir = cd;

% Normalise x2 and T to start at 0
if T(1) ~= 0 || x2(1) ~= 0
    disp('Data requires normalisation');
    x2n = x2 - x2(1);
    Txn = T - T(1);
else
    % Already normalised to remove transient
    x2n = x2;
    Txn = T;
    disp('Data already normalised');
end

% Obtain no. x2 births
nBirx2 = max(x2n);
% Clip number of x2 events to equal number of desired runs
if nBirx2  > runs
    id = find(x2n == runs, 1, 'last');
    x2n = x2n(1:id);
    % Set restricted Gillespie data set
    outComb.x2res = x2n;
    outComb.x1res = x1(1:id);
    outComb.Tres = Txn(1:id);
    outComb.restrict = 1;
elseif nBirx2 < runs
    % Correct the runs value if too many QBs specified and ensure that
    % later nFiles will be integer
    runsOrig = runs;
    runs = floor(nBirx2/runsdiv)*runsdiv;
    disp(['Runs modified from ' num2str(runsOrig) ' to ' num2str(runs)]);

    % Set restricted Gillespie data set
    id = find(x2n == runs, 1, 'last');
    x2n = x2n(1:id);
    outComb.x2res = x2n;
    outComb.x1res = x1(1:id);
    outComb.Tres = Txn(1:id);
    outComb.restrict = 1;
end

% Obtain reference times for photons which are set as the first time that
% the x2 number shows the increment (instead of last time at lower value,
% hence id = id + 1 == find([0 diff(x2n)] == 1)))
id = find(diff(x2n) == 1);
id = id + 1;
Tevent = Txn(id);
clear Txn x2n x1

% Check to ensure correct length of Tevent
if length(Tevent) ~= runs
    assignin('base', 'runs', runs);
    assignin('base', 'Tevent', Tevent);
    error('Tevent size and runs not consistent');
end

% Determine number of QB files to be loaded
nFiles = runs/runsdiv;
if nFiles ~= round(nFiles)
    error('Non-integer ratio of runs/runsdiv');
end

% Get time step for downsampling, assumes Nikolic data at same sample size
cd(folder.QBload);
% Expect QB files to be named as QBdiv[num]
load([fileRoot.loadQB '1'], 'tbio');
cd(thisDir);
% Assumes even sampling for tstep
tstep = tbio(2) - tbio(1);
timeStep = nStep*tstep;
disp(['Sampling with time steps of ' num2str(timeStep) ' ms']);

% Set interpolating start and end time for overall stream
Tstart = tbio(1) + Tevent(1);
Tend = Tevent(end) + tbio(length(tbio));
Tset = Tstart:timeStep:Tend;

% Declare variables for holding combined stream data
Icomp = zeros(size(Tset));
Ideclare = zeros(size(Tset));
Tqbdeclare = zeros(runsdiv, length(tbio));
Tstart = zeros(nFiles, runsdiv);
Tend = zeros(nFiles, runsdiv);

% Load QB files and shift their timeseries based on Tevent (photon times)
% and interpolate every QB current wave to this new time series with
% extrapolations set as 0
for i = 1:nFiles
    % Load a QB (divided) data file
    cd(folder.QBload);
    QBname1 = [fileRoot.loadQB num2str(i)];
    load(QBname1, 'Iset', 'tbio');
    cd(thisDir);

    % Declare some loop variables
    Iqb = Iset;
    clear('Iset');
    Tqb = Tqbdeclare;
    Iset = Ideclare;

    % Loop to obtain current stream after first offsetting the QBs with
    % photon times from the Gillespie simulation
    for j = 1:runsdiv
        % Declare and assign time values to account for event offsets
        Tqb(j, :) = tbio + Tevent(j + (i-1)*runsdiv);
        Tstart(i, j) = Tqb(j, 1);
        Tend(i, j) = Tqb(j, end);

        % Obtain interpolation range of times for each QB
        interpfirst = find(Tset >= Tstart(i, j), 1, 'first');
        interplast = find(Tset <= Tend(i, j), 1, 'last');
        interpRange = interpfirst:interplast;
        Tinterp = Tset(interpRange);

        % For each QB interpolate until last Tset in range then extrapolate
        % with zeros and sum into a stream, nakeinterp requires columns
        Iinterp = nakeinterp1(Tqb(j, :)', Iqb(j, :)', Tinterp');
        Icombine = Ideclare;
        Icombine(interpfirst:interplast) = Iinterp;
        Iset = Icombine + Iset;
    end

    % Save interpolated data files
    cd(folder.QBsave);
    QBname2 = [fileRoot.saveQB num2str(i)];
    save(QBname2, 'Tset', 'Iset');
    cd(thisDir);

    % Obtain complete overlapped stream
    Icomp = Icomp + Iset;
    clear('Iset');
    disp(['Completed interpolation for file' num2str(i)]);
end
% Completely combined stream after summing shifted bumps
Iset = Icomp;

% Assign some final outputs
outComb.runs = runs;
outComb.Tevent = Tevent;
outComb.Tstart = Tstart;
outComb.Tend = Tend;
outComb.Tset = Tset;
outComb.Iset = Iset;
outComb.timeStep = timeStep;