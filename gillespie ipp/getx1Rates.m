% Function that supplies rates based on the Q matrix of a Markov chain
% which represents x1's behaviour, basically a lookup
function rdotx1 = getx1Rates(Qx1, x1, nx1Reacs)

% Assumptions and modifications
% - accounts for non-nearest neighbour reactions
% - a birth death pair is obtained from Q for every reaction => r_const

% A lookup table on Q of x1 with the row corresponding to state x1 
% representing its reactions, x1 starts at 0 so increment by 1 for index
nState = length(Qx1);
state = x1 + 1;
if state < 1
    assignin('base', 'x1', x1);
    error('Incorrect state value');
end
Qrow = Qx1(state, 1:nState);

% Construct the rdot variable based on the Q matrix and current state by
% finding values left and right of the diagonal and accounting for
% dimensions of matrix in determining how many values exist
diagLoc = state;
rdotx1 = zeros(1, nx1Reacs);
if diagLoc + nx1Reacs/2 <= nState
    Qbir = Qrow(diagLoc+1:diagLoc+nx1Reacs/2);
else
    Qbir = Qrow(diagLoc+1:end);
end
if diagLoc - nx1Reacs/2 >= 1
    Qdea = Qrow(diagLoc-1:-1:diagLoc-nx1Reacs/2);
else
    Qdea = Qrow(diagLoc-1:-1:1);
end
Qbirlen = length(Qbir);
Qdealen = length(Qdea);

% Append zeros to account for all possible reactions
desiredLen = nx1Reacs/2;
if Qbirlen < desiredLen
    nPad = desiredLen - Qbirlen;
    Qbir = [Qbir zeros(1, nPad)];
end
if Qdealen < desiredLen 
    nPad = desiredLen - Qdealen;
    Qdea = [Qdea zeros(1, nPad)];
end

% Check consistent lengths and assign some check variables
if length(Qdea) ~= length(Qbir)
    error('The birth and death lengths are not consistent');
end

% Assign the rdot values with account for reactions which cannot occur at
% the given state and maintain [bir dea bir dea...] form
for i = 1:nx1Reacs
    if rem(i, 2) == 1
        % Odd reactions are births
        rdotx1(i) = Qbir((i+1)/2);
    else
        % Even reactions are deaths
        rdotx1(i) = Qdea(i/2);
    end
end