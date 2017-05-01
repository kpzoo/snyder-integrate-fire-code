function Iq = ghk_current(N,Cq_i,Cq_o,P1max,wq,z,V,S)

% Concentrations in mM
% P1max in um/sec
% Voltage in V
% Surface in um^2
% Current Iq in pA (>0 for in->out)

F = 96500;          % Faraday constant in C/mol
R = 8.31;           % Gas constant in J/(K*mol)
T = 300;            % Temperature in K

beta = z*F/(R*T);

Iq = 1e-6*N*wq*P1max*z*F*beta*V*S*( Cq_i - Cq_o*exp(-beta*V) )/( 1 - exp(-beta*V) );

