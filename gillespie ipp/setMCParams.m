% Sepcific function to generate general MCs to represent the x1 stochastic
% dynamics with the requirement that there is only one x2 reaction with
% form x2dot = alpha*x1 
function [molecType, crossType, reacType, bulk, transit] = setMCParams(jmax, r_const, MCform)

% Assumptions and modifications
% - only one type of MC form is available

% Obtain a pre-specified MC form
switch(MCform)
    case 1
        % The case where reactions of all sizes 1:jmax exist - can use a
        % zero r_const value to remove the effect of certain transitions
        lenr = length(r_const);
        
        % Obtain reaction and reactant types with all x1 reactions of linear
        % restricted space form
        molecType = [ones(1, lenr-1) 2];
        crossType = ones(1, lenr);
        reacType = [ones(1, lenr-1) 2];
        
        % Obtain reaction jump size based arrays
        bulk = [1:jmax 1:jmax];
        bulk = sort(bulk);
        transitx1 = [bulk 0];
        transitx2 = [zeros(1, lenr-1) 1];
        bulk = [bulk 1];
        transitx1(2:2:end) = -transitx1(2:2:end);
        transit = [transitx1; transitx2];         
        
    otherwise
        error(['Unsupported MC formulation for MCform = ' num2str(MCform)]);
end
        


