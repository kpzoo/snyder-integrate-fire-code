% Function to calculate MSE based on interpolated timeseries (the integrals
% use times from the ODE solver output)
function statSet = getStatsAltMSE(tcat, xcat, x, T)

% Assumptions and modifications
% - uses iterative methods to handle large datasets upon fine interpolation
% - xcat and tcat are outputs of Snyder ODEs and not evenly spaced
% - xcat is the estimated form of x at tcat, T is where x is defined
% - in calculating MSE uses evenly spaced samples at inter spacing

% Remove duplicate points - the +1 in id calculation chooses the latter of
% the duplicates, remove if want to remove the former duplicate
dtm = diff(tcat);
id = find(dtm == 0) + 1; % choose latter index
lent = length(tcat);

% Obtain altered data vectors with duplicate time indices removed
idset = 1:lent;
trunc = setdiff(idset', id);
tm = tcat(trunc);
xcapm = xcat(trunc);

% Get sample interval spacing and total sample points 
inter = mean(diff(tm))/10;
nSamps = (max(tm) - min(tm))/inter; % (non-integer)
% Get the no. samples in each set
nSets = 10000;
nSetSamps = floor(nSamps/nSets);
tSetSamps = (max(tm) - min(tm))/nSetSamps;

% Case in which both x1m and x1capm are to be linearly interpolated -
% perform iterative calculation of the statistics of interest
sum_em = 0;
sum_emSq = 0;
sum_len = 0;
tlim = -ones(1, nSetSamps+1);
tlim(1) = tm(1);  

% Loop across sets iteratively calculating statistics
for i = 2:nSetSamps
    % Define time limits and obtain ids of relevant section of data
    tlim(i) = tSetSamps + tlim(i-1);
    idlim = find(tm >= tlim(i-1) & tm < tlim(i));
    
    % Get sample times
    Tsamp = linspace(tm(idlim(1)), tm(idlim(end)), nSetSamps);
    % Truncate to relevant section
    Ttemp = tm(idlim);
    xcaptemp = xcapm(idlim);
    
    % Obtain sample times and interpolate 
    xSamp = interp1(T, x, Tsamp, 'previous');
    xcapSamp = interp1(Ttemp, xcaptemp, Tsamp);
    eSamp = xSamp - xcapSamp;
    
    lenSamp = length(eSamp);
    if lenSamp ~= length(Tsamp)
        assignin('base', 'Tsamp', Tsamp);
        assignin('base', 'eSamp', eSamp);
        error(['The sampling of the error curve failed at i = ' num2str(i)]);
    end
    
    % Obtain iterative sums for statistics - sampled evenly
    eSampSq = eSamp.*eSamp;
    sum_em = sum_em + trapz(Tsamp, eSamp);
    sum_emSq = sum_emSq + trapz(Tsamp, eSampSq);
end

% Check that the correct number of points were obtained and calculate means
if sum_len > nSamps 
    assignin('base', 'tlim', tlim);
    error(['Inconsistent sample size: [nSamps sum_len] = ' [num2str(nSamps)...
        ' ' num2str(sum_len)]]);
else
    % Get time averages of statistics
    em_mean = sum_em/range(tm);
    em_mse = sum_emSq/range(tm);
    em_var = em_mse - em_mean^2;
end

% Assign stats set with indicator string to indicate data order
statSet.name = {'mean', 'var', 'mse'};
statSet.vals = [em_mean em_var em_mse];