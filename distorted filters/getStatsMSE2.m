% Function that interpolates time-series for x1 and calculates MSE
function statSet = getStatsMSE2(Tset, xset, T, x)

% Assumptions and modifications
% - no Tevent needed 
% - assumes Tset has same value at Tset{i-1}(end) and Tset{i}(1)
% - samples will be uniformly spaced in time between events but not totally
% - there will have duplicate times due to jumps but unimportant for samples
% - true x is exact so zoh sampling

% Do sampling between event times, nSampEv each
nsampEv = 1000;
nEv = length(Tset);

% Variables to iteratively calculate statistics
sum_e = 0;
sum_eSq = 0;

% Loop between the sample times
for i = 1:nEv
    % Sample times of interest
    tsampEv = linspace(Tset{i}(1), Tset{i}(end), nsampEv);
    % True values, zoh sampling
    xsamp = interp1(T, x, tsampEv, 'previous');
    % Linearly interpolated samples from estimated function
    xhsamp = interp1(Tset{i}, xset{i}, tsampEv);
    % Error contribution to integrals
    esamp = xsamp - xhsamp;
    sum_e = sum_e + trapz(tsampEv, esamp);
    sum_eSq = sum_eSq + trapz(tsampEv, esamp.^2);
end

% Divide integrals by total time to get time average
dTev = Tset{end}(end) - Tset{1}(1);
e_mean = sum_e/dTev;
e_mse = sum_eSq/dTev;
e_var = e_mse - e_mean^2;

% Output calculated statistics
statSet.nsampEv = nsampEv;
statSet.nEv = nEv;
statSet.name = {'mean', 'var', 'mse'};
statSet.vals = [e_mean e_var e_mse];