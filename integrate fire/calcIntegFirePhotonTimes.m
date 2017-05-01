% Function to calculate the integrate-fire estimated photons - integrate a
% QB stream, use points where it crosses multiples of threshold => photons
function [Test, Iest, Qest, levels, Qinp] = calcIntegFirePhotonTimes(Iinp, Tinp, thresh)

% Assumptions and modifications
% - inputs of current (Iinp) - time (Tinp) for QBs and threshold (thresh)
% - outputs estimated times and current and integrated current at times

% Integrate current with respect to time
Qinp = cumtrapz(Tinp, Iinp);
Qmax = Qinp(end);

% Obtain nLevels as estimated no. photons 
nLevels = floor(Qmax/thresh);
levels = thresh*(1:nLevels);
% Use the threshold levels to extract photon time ids
idcross = zeros(1, nLevels);
for i = 1:nLevels
    % Index controlling each time a level is crossed
    idcross(i) = find(Qinp <= levels(i), 1, 'last');
end
% Extracted estimated location times, currents and current integrals
Test = Tinp(idcross);
Qest = Qinp(idcross);
Iest = Iinp(idcross);