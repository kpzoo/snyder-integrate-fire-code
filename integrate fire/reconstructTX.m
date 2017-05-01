% Reconstruct a new X vector from SSA data which has the new estimated
% photon times (x2 births) that are supplied
function [Tnew, Xnew] = reconstructTX(T, X, TphEst)

% Assumptions and modifications
% - remove repetitions in photon times by slight displacement

% Remove repetitions in estimates if they exist
[~, normVal, newlen, oldlen] = removeRepetition(TphEst);
if newlen == oldlen && normVal == 0
    % Case of no repetitions
    disp('No initial repetitions in photon stream');
else
    % Loop function until all repetitions removed
    while newlen ~= oldlen
        [TphEst , ~, newlen, oldlen] = removeRepetition(TphEst);
    end
end

% Use photon times to generate corresponding x2
x2new = 1:length(TphEst);
x2new = x2new';

% Obtain original x1 and x2
x1 = X(:, 1);
x2 = X(:, 2);

% Remove old x2 times from data vectors
idrem = find(diff([0; x2]) == 1);
idleft = setdiff(1:length(T), idrem);
Trem = T(idleft);
x2rem = x2(idleft);
x1rem = x1(idleft);

% Add new x2 times and sort T vector
Tnew = [Trem; TphEst];
[Tnew, idsort] = sort(Tnew);

% Apply sorting indices to x1 and x2 and reformat x2 and x1 suitably
x1new = [x1rem; zeros(size(Tnew))];
x1new = x1new(idsort);
x2new = [zeros(size(x2rem)); x2new];
x2new = x2new(idsort);
% assignin('base', 'x1newL', x1new);
% assignin('base', 'x2newL', x2new);

% Adjust values of x1new and x2new for consistency - x1new has zeros that
% need to be replaced with values  matching previous event and x1 has zeros
% that need to match surrounding values at the new photon times
idcurr = zeros(size(TphEst));
idcurr(1) = find(Tnew == TphEst(1));
for i = 2:length(TphEst) 
    % Obtain current photon event id
    idcurr(i) = find(Tnew == TphEst(i));
    
    % Indices before the id must have same x2 value as only at the event
    % id is x2 incremented
    if idcurr(i) ~= 1
        x2new(idcurr(i-1):idcurr(i)-1) = x2new(idcurr(i-1));
    end
    
    % Indices at the sort id must have previous value as no x1 event
    % occurred at this point (at i = 1 its zero)
    if idcurr(i-1) ~= 1
        x1new(idcurr(i-1)) = x1new(idcurr(i-1)-1);
    end
end

% Assign output
Xnew = [x1new x2new];
    

% % Test inputs to try on function
%  T = [0 1 3 7 15 16 17 20]';
%  x1 = [0 1 1 1 2 2 1 1]';
%  x2 = [0 0 1 2 2 3 3 4]';
%  X = [x1 x2];
%  TphEst = [10 100]';