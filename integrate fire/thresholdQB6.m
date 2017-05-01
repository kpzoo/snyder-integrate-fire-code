% Modification to allow extra input of the threshAvgRatio in order to
% determine if the spike metric truly matches the Snyder MSE

% Modification to include the threshold ratio in the test plots via extra
% input into the calcPhotEst subfunction
% Modification to output the ratio of thresh to avgLim in order to gauge
% how the optimal threshold links to the average QB integral

% Modification to the calcPhotonTimeEst function to allow for simple
% thresholding with threshMeth = 2
% Modification to subfunction calcPhotonTimeEst to produce test plots

% Modification to allow test metric to also be calculated across the
% complete threshold set simply for comparison with Snyder
% Modification to allow use of training and testing data

% Function that sequentially processes a QB current stream and learns an
% optimal threshold from training data before being used on test data to
% evaluate the accuracy of estimating photon times
function outThresh = thresholdQB6(Ttrain, Itrain, Ptrain, Ttest, Itest, Ptest,...
    threshMeth, rangeVal, rangeBool, dtcost)

% Assumptions and modifications

% Check that output data are photon times and not intensities
if ~(issorted(Ptrain) && issorted(Ptest))
    assignin('base', 'PtrainErr', Ptrain);
    assignin('base', 'PtestErr', Ptest);
    error('The output data are not an increasing sequence of photon times');
end

% Check that the training and test data are all 1D (no time histories are
% used here just a current-time stream to map to photon times)
if all(size(Itrain) > 1) || all(size(Itest) > 1) || all(size(Ttrain) > 1) || all(size(Ttest) > 1)
    error('The input is not simply a current-time stream');
end

% Normalise the current inputs to the range [0 1] and ensure normalisation
% is same for both training and test data
% assignin('base', 'ItrainOrig', Itrain);
% assignin('base', 'ItestOrig', Itest);
Itrain = -Itrain;
Itest = -Itest;
maxNorm = max(max(Itrain), max(Itest));
outThresh.maxNorm = maxNorm;
Itrain = Itrain/maxNorm;
Itest = Itest/maxNorm;
% assignin('base', 'Itrain', Itrain);
% assignin('base', 'Itest', Itest);
% assignin('base', 'Ttrain', Ttrain);
% assignin('base', 'Ttest', Ttest);
% assignin('base', 'Ptrain', Ptrain);
% assignin('base', 'Ptest', Ptest);

% Train via specified threshold detection method
if ~rangeBool
    % Default range for handling single or several threshold
    % methods with optimisation to a single metric value
    switch(threshMeth)
        case 1
            % Threshold gradient method
            %limVal = 0.0001:0.0005:0.02;
            limVal = 0.001:0.005:0.2;
        case 2
            % Pure thresholding
            limVal = 0.001:0.005:0.2;
        case 3
            % Integrate and fire - threshold calculated via average
            % integral multiplied by a factor
            Qmaxtemp = trapz(Ttrain, Itrain);
            avgLim = Qmaxtemp/length(Ptrain);
            limVal = avgLim*(0.2:0.1:2);
    end
else
    % Method runs a single threshold without optimisation - used
    % with a single threshold method
    limVal = rangeVal;
    if length(limVal) ~= 1 || length(threshMeth) ~= 1
        error('Incorrect limVal and threshMeth for the non-optimised method');
    end
end


% Declare variables and set plotting boolean
plotEst = 0;
dist = 0*limVal;
TtrainEst = cell(1, length(limVal));
TtrainReal = cell(1, length(limVal));

% Assign training inputs for use in loop optimisation
inpPh.Iinp = Itrain;
inpPh.Tinp = Ttrain;
inpPh.Pinp = Ptrain;
testPlot = 0;
    
% Loop through possible thresholds and determine the optimum value
% based on the spike metric
for i = 1:length(limVal)

    % Assign sub-function inputs in this training case
    inpPh.thresh = limVal(i);

    % Sub-function that calculates photon time estimates for a
    % given threshold value
    outPh = calcPhotonTimeEst(threshMeth, inpPh, testPlot, 0);
    TtrainEst{i} = outPh.Test;
    TtrainReal{i} = outPh.Treal;

    % Calculate the index of performance via a cost metric
%     dist(i) = spkd(TtrainEst{i}, TtrainReal{i}, dtcost);
    dist(i) = abs(length(TtrainEst{i}) - length(TtrainReal{i}));
end

% Use the minimum metric to determine threshold and assign output
[dmin imin] = min(dist);
thresh = limVal(imin);
threshAvgRatio = thresh/avgLim;
outThresh.distTrainSet = dist;
outThresh.distTrain = dmin;
outThresh.thresh = thresh;
outThresh.threshAvgRatio = threshAvgRatio;
assignin('base', 'avgLim', avgLim);

% Plot the best training estimate and the distance measure
if plotEst && ~rangeBool
    figure;
    plot(TtrainReal{imin}, 'bo');
    hold on
    plot(TtrainEst{imin}, 'ro');
    hold off
    xlabel('no. photons');
    ylabel('photon times');
    legend('real', 'estimate', 'location', 'best');
    title('Photon time estimates on training data');

    if threshMeth == 3
        figure;
        plot(limVal/avgLim, dist);
        xlabel('threshold value (normalised)');
        ylabel('metric distance');
        title('Metric behaviour across threshold range for normalised data');
    else
        figure;
        plot(limVal, dist);
        xlabel('threshold value');
        ylabel('metric distance');
        title('Metric behaviour across threshold range for normalised data');
    end
end

% Use threshold to estimate the photon times for the test data and
% assign to output structure
inpPh.Iinp = Itest;
inpPh.Tinp = Ttest;
inpPh.Pinp = Ptest;
inpPh.thresh = outThresh.thresh;
testPlot = 0;
outPh = calcPhotonTimeEst(threshMeth, inpPh, testPlot, threshAvgRatio);

% Calculate metric for test data and assign outputs
TtestReal = outPh.Treal;
TtestEst = outPh.Test;
% distTest = spkd(TtestEst, TtestReal, dtcost);
distTest = abs(length(TtestReal) - length(TtestEst));
outThresh.Treal = TtestReal;
outThresh.Test = TtestEst;
outThresh.distTest = distTest;

% Plot the performance on the test data
if plotEst && ~rangeBool
    figure;
    plot(TtestReal, 'bo');
    hold on
    plot(TtestEst, 'ro');
    hold off
    xlabel('no. photons');
    ylabel('photon times');
    legend('real', 'estimate', 'location', 'best');
    title('Photon time estimates on testing data');
end


% Sub-function to calculate photon times based on threshold, can be used on
% either training or testing data
function outPh = calcPhotonTimeEst(threshMeth, inpPh, testPlot, threshAvgRatio)

switch(threshMeth)
    case 1

        % Threshold based on signal level and gradient with a cost metric
        thresh = inpPh.thresh;
        Iinp = inpPh.Iinp;
        Tinp = inpPh.Tinp;
        Pinp = inpPh.Pinp;

        % Plot threshold method step by step if specified
        if testPlot
            figure;
            hold on
            plot(Tinp, Iinp);
            plot(Tinp, thresh*ones(size(Tinp)), 'k');
        end

        % Find points where current exceeds threshold
        id = find(Iinp > thresh);
        Iexcess = Iinp(id);
        Texcess = Tinp(id);
        if testPlot
            plot(Texcess, Iexcess, 'ro');
            assignin('base', 'Texcess', Texcess);
            assignin('base', 'Iexcess', Iexcess);
        end

        % Obtain successive indices and account for if idless < 1 or
        % idmore > length(Iinp)
        idless = id - 1;
        idmore = id + 1;
        irem = [find(idmore > length(Iinp)) find(idless < 1)];

        % Obtain difference in currents across the threshold points from
        % the successive indices
        if ~isempty(irem)
            % Remove the invalid entries for this analysis - should not
            % result in loss of photon as there will be nearby values
            % generally which also cross the threshold
            ileave = setdiff(1:length(id), irem);
            idless = idless(ileave);
            idmore = idmore(ileave);
            Iexcess = Iexcess(ileave);
            Texcess = Texcess(ileave);
        end
        Iless = Iinp(idless);
        Imore = Iinp(idmore);
        Idiff1 = Imore - Iexcess;
        Idiff2 = Iexcess - Iless;

        % Limit valid points under positive gradient requirement
        idlim = find(Idiff1 > 0 & Idiff2 > 0);
        Ilim = Iexcess(idlim);
        Tlim = Texcess(idlim);
        if testPlot
            plot(Tlim, Ilim, 'ms');
            assignin('base', 'Tlim', Tlim);
            assignin('base', 'Ilim', Ilim);
        end

        % Remove groups of sequential ids in idlim as these correspond
        % to the same QB i.e. all of its points which are above limVal
        idgr = find(diff(idlim) > 1);
        idgr = idgr + 1;
        Igr = Ilim(idgr);
        Tgr = Tlim(idgr);
        if testPlot
            plot(Tlim, Ilim, 'ms');
            assignin('base', 'Tgr', Tgr);
            assignin('base', 'Igr', Igr);
        end

        % Obtain estimated photon times and assign real photon times
        Test = Tgr;
        Treal = Pinp';
        outPh.Test = Test;
        outPh.Treal = Treal;

    case 2

        % Threshold based on signal level only with a cost metric and
        % represents the most primitive thresholding sigbnalling possible
        thresh = inpPh.thresh;
        Iinp = inpPh.Iinp;
        Tinp = inpPh.Tinp;
        Pinp = inpPh.Pinp;

        % Plot threshold method step by step if specified
        if testPlot
            figure;
            hold on
            plot(Tinp, Iinp);
            plot(Tinp, thresh*ones(size(Tinp)), 'k');
            %             assignin('base', 'Tinp', Tinp);
            %             assignin('base', 'Iinp', Iinp);
        end

        % Find points where current exceeds threshold
        id = find(Iinp > thresh);
        Iexcess = Iinp(id);
        Texcess = Tinp(id);
        if testPlot
            plot(Texcess, Iexcess, 'ro');
            assignin('base', 'Texcess', Texcess);
            assignin('base', 'Iexcess', Iexcess);
        end

        % Remove points with successive indices, method assumes that the QB
        % wave must dip below the threshold and then hit it again for a new
        % photon to be estimated
        idgr = find(diff(id) > 1);
        idgr = idgr + 1;
        Igr = Iexcess(idgr);
        Tgr = Texcess(idgr);
        if testPlot
            plot(Tgr, Igr, 'ms');
            assignin('base', 'Tgr', Tgr);
            assignin('base', 'Igr', Igr);
        end

        % Obtain estimated photon times and assign real photon times
        Test = Tgr;
        Treal = Pinp';
        outPh.Test = Test;
        outPh.Treal = Treal;

    case 3

        % Integrate and fire method with threshold which is applied to
        % integrated signal to obtain photon locations
        thresh = inpPh.thresh;
        Iinp = inpPh.Iinp;
        Tinp = inpPh.Tinp;
        Pinp = inpPh.Pinp;

        % Integrate current with respect to time
        Qinp = cumtrapz(Tinp, Iinp);
        Qmax = Qinp(end);
        if testPlot
            assignin('base', 'Tinp', Tinp);
            assignin('base', 'Qinp', Qinp);
            assignin('base', 'Iinp', Iinp);
        end

        % Obtain nLevels as estimated no, photons and use the threshold
        % multiples to obtain photon times
        nLevels = floor(Qmax/thresh);
        levels = thresh*(1:nLevels);
        idgr = zeros(1, nLevels);
        for i = 1:nLevels
            idgr(i) = find(Qinp <= levels(i), 1, 'last');
        end
        Tgr = Tinp(idgr);
        Qgr = Qinp(idgr);
        Igr = Iinp(idgr);
        if testPlot
            assignin('base', 'Tgr', Tgr);
            assignin('base', 'Qgr', Qgr);
            assignin('base', 'Igr', Igr);
        end

        % Obtain estimated photon times and assign real photon times
        Test = Tgr;
        Treal = Pinp';
        outPh.Test = Test;
        outPh.Treal = Treal;

        % Figure showing the threshold points on the current and integrated
        % current curves
        if testPlot
            figure;
            plot(Tinp, Iinp);
            hold on
            plot(Tgr, Igr, 'ko');
            hold off
            xlabel('time');
            ylabel('current');
            title(['Integrate-fire threshold points with ratio of ' num2str(threshAvgRatio)]);
            
            figure;
            plot(Tinp, Qinp);
            hold on
            plot(Tgr, Qgr, 'ko');
            hold off
            xlabel('time');
            ylabel('current');
            title(['Integrate-fire threshold points with ratio of ' num2str(threshAvgRatio)]);
        end
end