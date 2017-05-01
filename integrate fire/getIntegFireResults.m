% Function calculates the integrate-fire estimates using the avg QB area
function [Test, dist, thresh] = getIntegFireResults(Tset, Iset, Pset, testPlot)

% Assumptions and modifications
% - based on thresholdQB6
% - no training done here, use avg integral as threshold so avgRatio is 1

% Normalise the current inputs to the range [0 1]
Iset = -Iset;
Iset = Iset/max(Iset);

% Get average coded QB area and make it threshold, thresh
Qmax = trapz(Tset, Iset);
thresh = Qmax/length(Pset);
threshAvgRatio = 1;

% Get integrate-fire estimated photon time
[Test, Iest, Qest, levels, Qset] = calcIntegFirePhotonTimes(Iset, Tset, thresh);

% Calculate metric for estimate
dist = abs(length(Pset) - length(Test));

% Figure showing the threshold points for integrate-fire solution
if testPlot
    figure;
    subplot(2, 1, 1);
    plot(Tset, Iset, Test, Iest, 'ko');
    xlabel('time');
    ylabel('current');
    title(['Current threshold points, ratio = ' num2str(threshAvgRatio)]);
    subplot(2, 1, 2);
    plot(Tset, Qset, Test, Qest, 'ko');
    hold on
    h = gca;
    for i = 1:100:length(levels)
        plot(h.XLim, [levels(i) levels(i)], 'r');
    end
    hold off
    xlabel('time');
    ylabel('integrated current');
    title(['Integrated current threshold points, ratio = ' num2str(threshAvgRatio)]);
end

