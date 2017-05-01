% Function to re-sample the stochastic reaction timeseries
function [xn, tn] = getSamplesExtrapQuad(x, t, tn, extrapID)

% Assumptions and modifications
% - only quadtree type splitting
% - allows extrapolation with NaN or 0

% Declare xn and the bounds on t
xn = zeros(size(tn));
maxt = max(t);
mint = min(t);
len = length(tn);

% Separate t into quadrants and obtain the minimum limit of each
quad = floor(0.25*length(t));
tQ{1} = t(1:quad);
tQ{2} = t(quad+1:2*quad);
tQ{3} = t(2*quad+1:3*quad);
tQ{4} = t(3*quad+1:end);
minSet = [min(tQ{1}) min(tQ{2}) min(tQ{3}) min(tQ{4})];

% For every tn sample first find relevant quadrant and then find
% the lower nearest neighbour
for i = 1:len
    % Find relevant quadrant
    coord = sum(tn(i) >= minSet); % >= vs > so that coord = 0 not possible
    % Take value from last event as the function has not changed since
    id = find(tQ{coord} <= tn(i), 1, 'last');
    id = id + (coord - 1)*quad;
    xn(i) = x(id);
end

% Extrapolate the outside values
xn(tn < mint) = extrapID;
xn(tn > maxt) = extrapID;