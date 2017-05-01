% Modifications to reduce the execution time of getSamplesExtrap via data
% splitting via quadtree type searches

% Modification to allow extrapolation with chosen values such as NaN or 0
% as indicators of the extrapolation

% Simple function to re-sample the stochastic reaction timeseries - updated
% to do uniform sampling with a specified interval or to use the actual
% samples of time provided
function [xn tn] = getSamplesExtrap2(ver, samp_inter, x, t, tn, extrapID)

% Set method for executing interpolation
meth = 1;

% Obtain appropriate output vector length and calculate tn if needed
if ~ver
    % Divide t into number of samples if none provided
    tn = min(t):samp_inter:max(t);
    len = length(tn);
    disp(['Actual number of samples is ' num2str(len)]);
else
    len = length(tn);
end

% Declare xn and the bounds on t
xn = zeros(size(tn));
maxt = max(t);
mint = min(t);

switch(meth)
    case 0
        % Obtain relevant x samples noting that x is stepwise continuous
        for i = 1:len
            % Extrapolate if tn is outside range of t
            if tn(i) < mint || tn(i) > maxt
                xn(i) = extrapID;
            else
                % Take value from last event as the function has not changed since
                id = find(t <= tn(i), 1, 'last');
                xn(i) = x(id);
            end
        end
    case 1
        % Separate t into quadrants and obtain the minima
        quad = floor(0.25*length(t));
        tQ{1} = t(1:quad);
        tQ{2} = t(quad+1:2*quad);
        tQ{3} = t(2*quad+1:3*quad);
        tQ{4} = t(3*quad+1:end);
        minSet = [min(tQ{1}) min(tQ{2}) min(tQ{3}) min(tQ{4})];
        assignin('base', 'tQ', tQ);
        
        
        % For every tn sample first find relevant quadrant and then find
        % the lower nearest neighbour
        for i = 1:len
            % Find relevant quadrant
            coord = sum(tn(i) >= minSet); % <------- modified from just > to >= so that coord = 0 not possible
            % Take value from last event as the function has not changed since
            id = find(tQ{coord} <= tn(i), 1, 'last');
            id = id + (coord - 1)*quad;
            xn(i) = x(id);
        end
        
        % Extrapolate the outside values
        xn(tn < mint) = extrapID;
        xn(tn > maxt) = extrapID;
            
end