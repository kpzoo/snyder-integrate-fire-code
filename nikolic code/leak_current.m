function Iq = leak_current(z, D, Cmv, Cbody, L, S)

% Concentrations in mM
% Surface and length in um^2 and um
% D in um^2/sec
% Current Iq in pA

F = 96500;          % Faraday constant in C/mol

Iq = 1e-6*z*F*D*(S/L)*(Cmv - Cbody);

