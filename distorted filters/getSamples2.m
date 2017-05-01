% Simple function to re-sample the stochastic reaction timeseries - updated
% to do uniform sampling with a specified interval or to use the actual
% samples of time provided
function [xn tn] = getSamples2(ver, samp_inter, x, t, tn)

% Obtain appropriate output vector length and calculate tn if needed
if ~ver
    % Divide t into number of samples if none provided
    tn = min(t):samp_inter:max(t);
    len = length(tn);
    disp(['Actual number of samples is ' num2str(len)]);
else
    len = length(tn);
end

% Check dimensions of original inputs - assumes x is one dimensional
if ~all(size(x) == size(t))
    assignin('base', 'xFnInp', x);
    assignin('base', 'tFnInp', t);
    error('Function inputs have inconsistent dimensions');
end

% Obtain relevant sample x values noting that x is stepwise continuous
xn = zeros(1, len);
for i = 1:len
    id = find(t <= tn(i), 1, 'last');
    xn(i) = x(id);
end
    