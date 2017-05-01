% Check that the output rates are all members of the Q matrix entries
function check = checkRateQ(Q, Xdot)

% Assumptions and modifications
% - only the last reaction is for x2 (standard throughout code)
% - assumes zeros for non-reactions in rdot (e.g if jumps of 1 and 3)

% Extract the x1 reactions and initialise check variable
rdot = Xdot(:, 1:end-1);
nReacs = size(rdot, 2);
check = zeros(1, nReacs);

% Noting the structure [bir dea bir dea]
for i = 1:nReacs
    
    % Obtain correct diagonal of Q matrix
    rdiv = rem(i, 2);
    if rdiv == 1
        % Odd reactions are births
        reac = (i+1)/2;
        Qdiag = diag(Q, reac);      
    else
        % Even reactions are deaths
        reac = i/2;
        Qdiag = diag(Q, -reac);
    end
    
    % Check the Q entries with the rate values (include 0 rate which is not
    % in the Q matrix)
    rateVal = unique(rdot(:, i));
    if ismember(rateVal, [Qdiag; 0])
        check(i) = 1;
    end
end

% Determine overall correctness
sumCheck = sum(check);
if sumCheck == nReacs
    disp('All reactions have the correct values from the Q matrix');
else
    disp(['There are ' num2str(nReacs - sumCheck) ' erroneous reactions']);
end
    
    